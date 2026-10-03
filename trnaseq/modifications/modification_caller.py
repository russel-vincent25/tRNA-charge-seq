"""
Modification Caller for tRNA RT Signatures

This module annotates RT signatures with specific tRNA modification types based on
known RT signature patterns from the literature.

Sites are detected by OR across three channels -- substitution, deletion and
RT stop -- each tested against its own beta-binomial background, because the
channel a modification surfaces in depends on the reverse transcriptase.
Modification identity is assigned afterwards as a ranked candidate list
(see docs/CHANNEL_COMBINED_DESIGN.md).

Key modifications detected:
- m1A (1-methyladenosine): Strong RT stops, A->any mismatches at positions 58, 14
- m3C (3-methylcytosine): C->T mismatches at positions 32
- Ψ (pseudouridine): U->C signature at multiple positions
- m7G (7-methylguanosine): G->C/A signature at position 46
- m5C (5-methylcytosine): Subtle C->T signature at positions 48, 49
- i6A (N6-isopentenyladenosine): A->G at position 37

References:
- Carlile et al. 2014 (Nature) - Pseudouridine detection
- Schwartz et al. 2014 (Cell) - m1A, m3C detection
- Hauenschild et al. 2015 (NAR) - Modification signatures
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, replace
from scipy import stats
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.stats import binomtest, combine_pvalues

from .channel_priors import CHANNELS, DERIVATION_THRESHOLDS, load_channel_priors

SIGNATURE_TYPES = ('mismatch', 'rt_stop', 'gap', 'combined')

# Detection channels each signature_type may be called from
_SIGNATURE_CHANNELS = {
    'mismatch': ('mismatch',),
    'gap': ('deletion',),
    'rt_stop': ('rt_stop',),
    'combined': CHANNELS,
}

# Per-channel rate thresholds for calling (RT stop as a fraction). Same
# values the channel priors were derived at.
DEFAULT_CHANNEL_THRESHOLDS = dict(DERIVATION_THRESHOLDS)

# Signature-table column and scale giving each channel's rate
_RATE_SOURCE = {
    'mismatch': ('mismatch_rate', 1.0),
    'deletion': ('gap_rate', 1.0),
    'rt_stop': ('rt_stop_pct', 0.01),
}

# A fired channel whose prior weight is below this is 'unexpected' for the enzyme
UNEXPECTED_CHANNEL_WEIGHT = 0.10

_NT_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


@dataclass
class ModificationProfile:
    """
    Profile for a specific tRNA modification type.

    Attributes:
        name: Modification name (e.g., 'm1A')
        full_name: Full chemical name
        typical_positions: Common positions where this modification occurs
        signature_type: Type of RT signature ('mismatch', 'rt_stop', 'gap',
            'combined'). Defaults to 'combined' because the channel a
            modification surfaces in is set by the reverse transcriptase,
            not the modification; override per profile (or via
            ``ModificationCaller(signature_overrides=...)``) only when a
            channel is known to be uninformative.
        mismatch_pattern: Expected mismatch pattern (e.g., 'A->G', 'A->any')
        min_rate: Legacy per-profile mismatch floor; detection now uses the
            per-channel thresholds in ModificationCaller (unused for calling)
        rt_stop_required: Whether RT stops are required for calling
        min_rt_stop_pct: Minimum RT stop percentage if required
    """
    name: str
    full_name: str
    typical_positions: List[int]
    signature_type: str = 'combined'
    mismatch_pattern: Optional[str] = None
    min_rate: float = 0.10
    rt_stop_required: bool = False
    min_rt_stop_pct: float = 15.0
    confidence_weight: float = 1.0


# Known modification profiles based on literature
MODIFICATION_PROFILES = {
    'm1A': ModificationProfile(
        name='m1A',
        full_name='1-methyladenosine',
        typical_positions=[58, 14, 9],
        mismatch_pattern='A->any',
        min_rate=0.10,
        rt_stop_required=False,  # RT stop boosts confidence but isn't required
        min_rt_stop_pct=15.0,
        confidence_weight=1.5  # High confidence when both mismatch + RT stop
    ),

    'm3C': ModificationProfile(
        name='m3C',
        full_name='3-methylcytosine',
        typical_positions=[32],
        mismatch_pattern='C->T',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.2
    ),

    'pseudouridine': ModificationProfile(
        name='Ψ',
        full_name='pseudouridine',
        typical_positions=[27, 28, 31, 32, 39, 40, 55, 13, 38],
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm7G': ModificationProfile(
        name='m7G',
        full_name='7-methylguanosine',
        typical_positions=[46],
        mismatch_pattern='G->any',  # TGIRT/Maxima produce G->C; some enzymes G->A
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.1
    ),

    'm5C': ModificationProfile(
        name='m5C',
        full_name='5-methylcytosine',
        typical_positions=[48, 49, 34, 40],
        mismatch_pattern='C->T',
        min_rate=0.05,  # Subtle signature
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'i6A': ModificationProfile(
        name='i6A',
        full_name='N6-isopentenyladenosine',
        typical_positions=[37],
        mismatch_pattern='A->G',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm2G': ModificationProfile(
        name='m2G',
        full_name='N2-methylguanosine',
        typical_positions=[10, 26],
        mismatch_pattern='G->A',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    'm22G': ModificationProfile(
        name='m22G',
        full_name='N2,N2-dimethylguanosine',
        typical_positions=[26],
        mismatch_pattern='G->A',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    's4U': ModificationProfile(
        name='s4U',
        full_name='4-thiouridine',
        typical_positions=[8, 9, 4],
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm1G': ModificationProfile(
        name='m1G',
        full_name='1-methylguanosine',
        typical_positions=[37, 9],
        mismatch_pattern='G->any',
        min_rate=0.08,
        rt_stop_required=False,
        min_rt_stop_pct=15.0,
        confidence_weight=1.1
    ),

    'cmo5U': ModificationProfile(
        name='cmo5U',
        full_name='uridine 5-oxyacetic acid',
        typical_positions=[34],
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    'mnm5s2U': ModificationProfile(
        name='mnm5s2U',
        full_name='5-methylaminomethyl-2-thiouridine',
        typical_positions=[34],
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    't6A': ModificationProfile(
        name='t6A',
        full_name='N6-threonylcarbamoyladenosine',
        typical_positions=[37],
        mismatch_pattern='A->T',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'ms2i6A': ModificationProfile(
        name='ms2i6A',
        full_name='2-methylthio-N6-isopentenyladenosine',
        typical_positions=[37],
        mismatch_pattern='A->G',
        min_rate=0.08,
        rt_stop_required=False,
        min_rt_stop_pct=10.0,
        confidence_weight=1.1
    ),

    'I': ModificationProfile(
        name='I',
        full_name='inosine',
        typical_positions=[34],
        mismatch_pattern='A->G',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    # --- Eukaryotic-enriched modifications ---

    'ac4C': ModificationProfile(
        name='ac4C',
        full_name='N4-acetylcytidine',
        typical_positions=[12, 34],
        mismatch_pattern='C->T',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    'Gm': ModificationProfile(
        name='Gm',
        full_name="2'-O-methylguanosine",
        typical_positions=[18, 34],
        mismatch_pattern='G->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'Cm': ModificationProfile(
        name='Cm',
        full_name="2'-O-methylcytidine",
        typical_positions=[32, 34],
        mismatch_pattern='C->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'Um': ModificationProfile(
        name='Um',
        full_name="2'-O-methyluridine",
        typical_positions=[32, 44],
        mismatch_pattern='U->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'Am': ModificationProfile(
        name='Am',
        full_name="2'-O-methyladenosine",
        typical_positions=[44],
        mismatch_pattern='A->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'Q': ModificationProfile(
        name='Q',
        full_name='queuosine',
        typical_positions=[34],
        mismatch_pattern='G->any',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm6A': ModificationProfile(
        name='m6A',
        full_name='N6-methyladenosine',
        typical_positions=[37, 58],
        mismatch_pattern='A->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.9
    ),
}


def estimate_background_error_rate(
    pscm_dict: Dict[str, np.ndarray],
    ref_dict: Dict[str, dict],
    synthetic_prefixes: Tuple[str, ...] = ('Synthetic_',),
    min_coverage: int = 50,
) -> Tuple[float, str]:
    """Estimate background sequencing error rate from PSCM data.

    Strategy:
    1. If synthetic (spike-in) tRNAs are present, use their mismatch rates
       (weighted by coverage) as the background — these have no modifications.
    2. Fallback: compute per-position mismatch rate across ALL tRNAs and
       take the 25th percentile (most positions are unmodified).
    3. Floor the result at 0.001 to avoid zero-division in fold-change.

    Args:
        pscm_dict: {trna_name: ndarray(ref_len, 8)} — columns are
            A, C, G, T, N, gap, coverage, rt_stop.
        ref_dict: {trna_name: {'seq': str, 'seq_len': int}}.
        synthetic_prefixes: FASTA name prefixes that identify synthetic tRNAs.
        min_coverage: Ignore positions with coverage below this threshold.

    Returns:
        (error_rate, source) where source is 'synthetic' or 'empirical_q25'.
    """
    _NT_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    FLOOR = 0.001

    # --- Try synthetic tRNAs first ---
    total_mismatches = 0
    total_coverage = 0
    for trna_name, mat in pscm_dict.items():
        if not any(trna_name.startswith(p) for p in synthetic_prefixes):
            continue
        if trna_name not in ref_dict:
            continue
        ref_seq = ref_dict[trna_name]['seq'].upper()
        for pos_idx in range(mat.shape[0]):
            cov = mat[pos_idx, 6]  # coverage column
            if cov < min_coverage:
                continue
            ref_nt = ref_seq[pos_idx] if pos_idx < len(ref_seq) else 'N'
            if ref_nt not in _NT_IDX:
                continue
            correct = mat[pos_idx, _NT_IDX[ref_nt]]
            mis = cov - correct
            total_mismatches += max(0, mis)
            total_coverage += cov

    if total_coverage > 0:
        rate = total_mismatches / total_coverage
        return (max(FLOOR, rate), 'synthetic')

    # --- Fallback: 25th percentile of all positions ---
    rates = []
    for trna_name, mat in pscm_dict.items():
        if trna_name not in ref_dict:
            continue
        ref_seq = ref_dict[trna_name]['seq'].upper()
        for pos_idx in range(mat.shape[0]):
            cov = mat[pos_idx, 6]
            if cov < min_coverage:
                continue
            ref_nt = ref_seq[pos_idx] if pos_idx < len(ref_seq) else 'N'
            if ref_nt not in _NT_IDX:
                continue
            correct = mat[pos_idx, _NT_IDX[ref_nt]]
            mis = cov - correct
            rates.append(max(0, mis) / cov)

    if rates:
        q25 = float(np.percentile(rates, 25))
        return (max(FLOOR, q25), 'empirical_q25')

    return (FLOOR, 'empirical_q25')


# ---------------------------------------------------------------------------
# Per-channel background model
# ---------------------------------------------------------------------------

@dataclass
class ChannelBackground:
    """Null model for one detection channel.

    Beta-binomial with mean *mean* and overdispersion (intra-class
    correlation) *rho*; ``rho == 0`` is a plain binomial. A single binomial
    rate is not calibrated for any channel -- per-position null rates are
    overdispersed, RT stops ~10x more than mismatches -- so *rho* matters.

    Attributes:
        mean: Null event rate per read (substitution, deletion or RT stop).
        rho: Overdispersion in [0, 1); 0 means binomial.
        source: 'synthetic', 'empirical_bulk' or 'fixed'.
        n_positions: Positions the fit used.
    """
    mean: float
    rho: float = 0.0
    source: str = 'fixed'
    n_positions: int = 0

    def sf(self, k, n) -> np.ndarray:
        """P(X >= k) under the null, vectorised over positions."""
        k = np.asarray(k, dtype=np.int64)
        n = np.asarray(n, dtype=np.int64)
        if self.mean <= 0:
            return np.where(k > 0, 0.0, 1.0)
        if self.rho <= 1e-9:
            return stats.binom.sf(k - 1, n, self.mean)
        a = self.mean * (1 - self.rho) / self.rho
        b = (1 - self.mean) * (1 - self.rho) / self.rho
        return stats.betabinom.sf(k - 1, n, a, b)


def _channel_counts(mat: np.ndarray, ref_seq: str):
    """Per-position (k, n, valid) for each channel from a PositionalExtractor
    matrix (columns A, C, G, T, N, gap, coverage, rt_stop)."""
    n_pos = min(mat.shape[0], len(ref_seq))
    mat = mat[:n_pos]
    nt_idx = np.array([_NT_IDX.get(nt, -1) for nt in ref_seq.upper()[:n_pos]])
    valid = nt_idx >= 0
    cov = mat[:, 6]
    correct = np.where(valid, mat[np.arange(n_pos), np.clip(nt_idx, 0, 3)], 0)
    gap = mat[:, 5]
    counts = {
        'mismatch': np.clip(cov - correct - gap, 0, None),
        'deletion': gap,
        'rt_stop': mat[:, 7],
    }
    # Position 1: every full-length read "stops" there by construction
    stop_valid = valid.copy()
    stop_valid[:1] = False
    masks = {'mismatch': valid, 'deletion': valid, 'rt_stop': stop_valid}
    return cov, counts, masks


def _fit_beta_binomial(k: np.ndarray, n: np.ndarray) -> Tuple[float, float]:
    """ML fit of (mean, rho) for a beta-binomial; falls back to binomial."""
    total_k, total_n = float(k.sum()), float(n.sum())
    if total_k == 0:
        # No events at all: half a pseudo-event keeps p-values finite
        return 0.5 / total_n, 0.0
    mu0 = total_k / total_n

    def nll(t):
        mu, rho = expit(t[0]), expit(t[1])
        a = mu * (1 - rho) / rho
        b = (1 - mu) * (1 - rho) / rho
        return -stats.betabinom.logpmf(k, n, a, b).sum()

    res = minimize(nll, [logit(mu0), logit(1e-3)], method='Nelder-Mead',
                   options={'xatol': 1e-4, 'fatol': 1e-3, 'maxiter': 400})
    if not np.isfinite(res.fun):
        return mu0, 0.0
    return float(expit(res.x[0])), float(expit(res.x[1]))


def estimate_channel_backgrounds(
    pscm_dict: Dict[str, np.ndarray],
    ref_dict: Dict[str, dict],
    synthetic_prefixes: Tuple[str, ...] = ('Synthetic_',),
    min_coverage: int = 50,
    thresholds: Optional[Dict[str, float]] = None,
    max_positions: int = 5000,
    seed: int = 0,
) -> Dict[str, ChannelBackground]:
    """Fit a beta-binomial null per detection channel for one sample.

    Positions come from synthetic (unmodified) tRNAs when present, otherwise
    from all tRNAs. Either way only the null bulk is used -- positions whose
    rate is below half the channel's calling threshold -- so real
    modifications do not inflate the null. At most *max_positions* positions
    per channel are used (random, seeded) to bound fitting time.

    Args:
        pscm_dict: {trna_name: ndarray(ref_len, 8)} for ONE sample.
        ref_dict: {trna_name: {'seq': str, ...}}.
        synthetic_prefixes: Name prefixes of synthetic spike-in tRNAs.
        min_coverage: Ignore positions below this coverage.
        thresholds: Per-channel calling thresholds (default
            :data:`DEFAULT_CHANNEL_THRESHOLDS`).
        max_positions: Cap on positions per channel fit.
        seed: RNG seed for the subsample.

    Returns:
        {channel: ChannelBackground}.
    """
    thresholds = thresholds or DEFAULT_CHANNEL_THRESHOLDS

    def collect(names):
        ks = {c: [] for c in CHANNELS}
        ns = {c: [] for c in CHANNELS}
        for name in names:
            if name not in ref_dict:
                continue
            cov, counts, masks = _channel_counts(pscm_dict[name], ref_dict[name]['seq'])
            for c in CHANNELS:
                keep = masks[c] & (cov >= min_coverage)
                ks[c].append(counts[c][keep])
                ns[c].append(cov[keep])
        return ({c: np.concatenate(ks[c]) if ks[c] else np.array([]) for c in CHANNELS},
                {c: np.concatenate(ns[c]) if ns[c] else np.array([]) for c in CHANNELS})

    synthetic = [t for t in pscm_dict if t.startswith(tuple(synthetic_prefixes))]
    k_all, n_all = collect(synthetic)
    source = 'synthetic'
    if not any(len(n_all[c]) for c in CHANNELS):
        k_all, n_all = collect(list(pscm_dict))
        source = 'empirical_bulk'

    rng = np.random.default_rng(seed)
    backgrounds = {}
    for c in CHANNELS:
        k, n = k_all[c], n_all[c]
        if len(n) == 0:
            backgrounds[c] = ChannelBackground(mean=0.001, source='default')
            continue
        bulk = (k / n) < 0.5 * thresholds[c]
        k, n = k[bulk].astype(np.int64), n[bulk].astype(np.int64)
        if len(n) > max_positions:
            pick = rng.choice(len(n), max_positions, replace=False)
            k, n = k[pick], n[pick]
        mean, rho = _fit_beta_binomial(k, n)
        backgrounds[c] = ChannelBackground(mean=mean, rho=rho, source=source,
                                           n_positions=len(n))
    return backgrounds


def _bh_qvalues(pvals: np.ndarray, n_tests: int) -> np.ndarray:
    """BH q-values for *pvals* within a family of *n_tests* tests.

    Tests not passed in (sites below every effect threshold) are treated as
    having larger p-values than any passed in, which is conservative.
    """
    p = np.asarray(pvals, dtype=np.float64)
    m = max(int(n_tests), len(p))
    if len(p) == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.minimum(ranked, 1.0)
    return q


def _same_nt(a: str, b: str) -> bool:
    return a == b or {a, b} <= {'U', 'T'}


class ModificationCaller:
    """
    Call tRNA modification sites from per-position RT signatures.

    **Detection** is per site and channel-agnostic: a site is called when
    any of the three channels -- substitution (``mismatch``), deletion and
    RT stop -- *fires*. A channel fires when its rate is at or above that
    channel's threshold **and** (with ``statistical_test``) its p-value
    against that channel's own background survives Benjamini-Hochberg over
    all site x channel tests. Neither the modification profiles, the
    position priors nor the RT-enzyme channel priors influence detection.

    **Identity** is a separate step: every profile compatible with the
    site (reference nucleotide, substitution pattern when the mismatch
    channel fired, ``signature_type`` channel restriction) becomes a
    candidate, ranked by typical-position match, pattern fraction and
    ``confidence_weight``. One row per site; the top candidate is
    ``modification`` and the full ranked list is ``candidates``.

    **Confidence** is a ranking score, not a gate (``min_confidence``
    defaults to 0). Per-enzyme channel priors scale it, and flag
    ``unexpected_channel`` when signal arrives only in channels the enzyme
    rarely uses -- the call is kept either way.

    Example:
        >>> caller = ModificationCaller(organism='human', rt_enzyme='Maxima',
        ...                             channel_backgrounds=backgrounds)
        >>> calls = caller.call_all(trna_name, analysis['signatures'],
        ...                         pscm_df, ref_seq)
    """

    def __init__(
        self,
        organism: str = 'human',
        min_confidence: float = 0.0,
        use_position_priors: bool = True,
        statistical_test: bool = True,
        alpha: float = 0.01,
        background_error_rate: float = 0.01,
        rt_enzyme: Optional[str] = None,
        rt_temp: Optional[float] = None,
        signature_overrides: Optional[Dict[str, str]] = None,
        channel_backgrounds: Optional[Dict[str, ChannelBackground]] = None,
        channel_thresholds: Optional[Dict[str, float]] = None,
    ):
        """
        Initialize modification caller.

        Args:
            organism: Organism name (for position priors)
            min_confidence: Drop calls whose confidence (before channel-prior
                scaling) is below this. Default 0: detection is decided by
                the channel tests, not by confidence.
            use_position_priors: Use known modification positions to rank
                candidates and boost confidence
            statistical_test: Require per-channel significance (BH) in
                addition to the rate threshold
            alpha: FDR level for the per-channel tests
            background_error_rate: Legacy single null rate, applied to every
                channel as a binomial when *channel_backgrounds* is None.
            rt_enzyme: Reverse transcriptase that produced the data (e.g.
                'Maxima', 'SSIV', 'TGIRT', 'Indura'). Recorded on every call
                for provenance and used to look up per-enzyme channel
                priors. ``None`` means unknown: enzyme-agnostic behaviour.
            rt_temp: RT incubation temperature in degrees C, or ``None``.
            signature_overrides: ``{profile_key: signature_type}`` to
                restrict specific modifications to one channel (keys as in
                :data:`MODIFICATION_PROFILES`).
            channel_backgrounds: ``{channel: ChannelBackground}`` from
                :func:`estimate_channel_backgrounds`.
            channel_thresholds: Per-channel rate thresholds (RT stop as a
                fraction); default :data:`DEFAULT_CHANNEL_THRESHOLDS`.
        """
        self.organism = organism
        self.min_confidence = min_confidence
        self.use_position_priors = use_position_priors
        self.statistical_test = statistical_test
        self.alpha = alpha
        self.background_error_rate = background_error_rate
        self.rt_enzyme = rt_enzyme
        self.rt_temp = rt_temp
        self.channel_priors = load_channel_priors(rt_enzyme, rt_temp)

        if channel_backgrounds is None:
            channel_backgrounds = {
                c: ChannelBackground(mean=background_error_rate)
                for c in CHANNELS
            }
        missing = set(CHANNELS) - set(channel_backgrounds)
        if missing:
            raise ValueError(f"channel_backgrounds missing {sorted(missing)}")
        self.channel_backgrounds = channel_backgrounds
        self.channel_thresholds = dict(DEFAULT_CHANNEL_THRESHOLDS)
        self.channel_thresholds.update(channel_thresholds or {})

        # Load modification profiles (copied when overridden so the
        # module-level defaults are never mutated)
        self.profiles = dict(MODIFICATION_PROFILES)
        for key, sig_type in (signature_overrides or {}).items():
            if key not in self.profiles:
                raise ValueError(
                    f"Unknown modification profile '{key}' in "
                    f"signature_overrides; known: {sorted(self.profiles)}"
                )
            if sig_type not in SIGNATURE_TYPES:
                raise ValueError(
                    f"Invalid signature_type '{sig_type}' for '{key}'; "
                    f"expected one of {SIGNATURE_TYPES}"
                )
            self.profiles[key] = replace(self.profiles[key],
                                         signature_type=sig_type)

    def match_mismatch_pattern(
        self,
        pattern: str,
        ref_nt: str,
        pscm_row: pd.Series
    ) -> Tuple[bool, float]:
        """
        Check if observed substitution pattern matches expected pattern.

        Deletions ('-') and N calls are not substitutions and are excluded
        from the mismatch total, so a deletion-dominated site cannot fail
        or pass a substitution pattern on the strength of its gaps.

        Args:
            pattern: Expected pattern (e.g., 'A->G', 'A->any', 'C->T')
            ref_nt: Reference nucleotide
            pscm_row: Row from PSCM with nucleotide counts

        Returns:
            Tuple of (matches, fraction_of_pattern)
        """
        if pattern is None:
            return True, 0.0

        parts = pattern.split('->')
        if len(parts) != 2:
            return False, 0.0

        expected_ref, expected_obs = parts

        # Check reference nucleotide matches (handle U/T equivalence)
        if not _same_nt(expected_ref, ref_nt):
            return False, 0.0

        total_coverage = pscm_row.sum()
        if total_coverage == 0:
            return False, 0.0

        correct_count = pscm_row.get(ref_nt, 0)
        non_substitution = pscm_row.get('-', 0) + pscm_row.get('N', 0)
        total_mismatches = total_coverage - correct_count - non_substitution

        if total_mismatches <= 0:
            return False, 0.0

        # Check if specific nucleotide pattern matches
        if expected_obs == 'any':
            # Any substitution is acceptable
            fraction = total_mismatches / total_coverage
            return True, fraction
        else:
            # Specific nucleotide required (handle U/T equivalence)
            obs_count = pscm_row.get(expected_obs, 0)
            if obs_count == 0 and expected_obs in ('U', 'T'):
                alt = 'T' if expected_obs == 'U' else 'U'
                obs_count = pscm_row.get(alt, 0)
            if obs_count == 0:
                return False, 0.0
            fraction = obs_count / total_coverage
            # Pattern matches if this specific substitution is dominant
            matches = obs_count >= (total_mismatches * 0.5)
            return matches, fraction

    def perform_statistical_test(
        self,
        coverage: int,
        mismatch_count: int,
        expected_error_rate: float = None,
    ) -> float:
        """
        Binomial test of a single count against a fixed rate.

        Legacy helper; site calling uses the per-channel
        :class:`ChannelBackground` models instead.

        Args:
            coverage: Total read coverage
            mismatch_count: Number of events observed
            expected_error_rate: Expected rate. If *None*, uses
                ``self.background_error_rate``.

        Returns:
            P-value from binomial test
        """
        if expected_error_rate is None:
            expected_error_rate = self.background_error_rate
        if coverage == 0:
            return 1.0

        result = binomtest(
            k=int(mismatch_count),
            n=int(coverage),
            p=expected_error_rate,
            alternative='greater'
        )

        return result.pvalue

    # ------------------------------------------------------------------
    # Detection
    # ------------------------------------------------------------------

    @staticmethod
    def _tested_mask(signatures_df: pd.DataFrame, min_coverage: int) -> np.ndarray:
        if signatures_df is None or signatures_df.empty or 'coverage' not in signatures_df:
            return np.zeros(0, dtype=bool)
        cov = pd.to_numeric(signatures_df['coverage'], errors='coerce').fillna(0)
        return (cov >= max(1, min_coverage)).to_numpy()

    def count_tests(self, signatures_df: pd.DataFrame, min_coverage: int = 50) -> int:
        """Number of site x channel tests :meth:`call_all` performs.

        Sum this over every tRNA of a sample and pass it to
        :meth:`finalize_calls` for sample-wide FDR control.
        """
        tested = self._tested_mask(signatures_df, min_coverage)
        if not tested.any():
            return 0
        pos1 = int((signatures_df['position'].to_numpy()[tested] == 1).sum())
        return len(CHANNELS) * int(tested.sum()) - pos1

    def _test_channels(self, signatures_df: pd.DataFrame, min_coverage: int) -> pd.DataFrame:
        """Per-channel rate, effect gate and p-value at every tested site."""
        tested = self._tested_mask(signatures_df, min_coverage)
        if not tested.any():
            return pd.DataFrame()
        sub = signatures_df.loc[tested]
        positions = sub['position'].astype(int).to_numpy()
        coverage = pd.to_numeric(sub['coverage']).to_numpy(dtype=np.float64)
        n = np.rint(coverage).astype(np.int64)

        out = pd.DataFrame({'position': positions, 'coverage': coverage})
        for c in CHANNELS:
            col, scale = _RATE_SOURCE[c]
            rate = (pd.to_numeric(sub[col], errors='coerce').fillna(0).to_numpy(dtype=np.float64)
                    if col in sub.columns else np.zeros(len(sub))) * scale
            if c == 'rt_stop':
                rate = np.where(positions == 1, 0.0, rate)
            out[f'rate_{c}'] = rate
            out[f'gate_{c}'] = rate >= self.channel_thresholds[c]
            if self.statistical_test:
                k = np.rint(rate * n).astype(np.int64)
                p = self.channel_backgrounds[c].sf(k, n)
                if c == 'rt_stop':
                    p = np.where(positions == 1, 1.0, p)
                out[f'pvalue_{c}'] = p
            else:
                out[f'pvalue_{c}'] = np.nan
        return out

    # ------------------------------------------------------------------
    # Identity
    # ------------------------------------------------------------------

    def _rank_candidates(
        self,
        position: int,
        ref_nt: Optional[str],
        pscm_row: Optional[pd.Series],
        gated: Dict[str, bool],
        rt_stop_rate: float,
    ) -> List[Tuple[ModificationProfile, float, bool]]:
        """Compatible profiles at a site, best first: (profile, pattern_fraction, typical)."""
        ranked = []
        for profile in self.profiles.values():
            if not any(gated[c] for c in _SIGNATURE_CHANNELS[profile.signature_type]):
                continue
            if profile.rt_stop_required and rt_stop_rate * 100 < profile.min_rt_stop_pct:
                continue
            fraction = 0.0
            if profile.mismatch_pattern and ref_nt is not None:
                expected_ref = profile.mismatch_pattern.split('->')[0]
                if not _same_nt(expected_ref, ref_nt):
                    continue
                # Substitution identity is only observable when that channel fired
                if gated['mismatch'] and pscm_row is not None:
                    matches, fraction = self.match_mismatch_pattern(
                        profile.mismatch_pattern, ref_nt, pscm_row)
                    if not matches:
                        continue
            typical = self.use_position_priors and position in profile.typical_positions
            ranked.append((profile, fraction, typical))
        ranked.sort(key=lambda t: (t[2], t[1], t[0].confidence_weight), reverse=True)
        return ranked

    @staticmethod
    def _dominant_substitution(ref_nt: Optional[str], pscm_row: Optional[pd.Series]) -> str:
        if pscm_row is None or ref_nt is None:
            return ''
        counts = {f'{ref_nt}->{nt}': pscm_row.get(nt, 0)
                  for nt in ('A', 'C', 'G', 'T') if nt != ref_nt}
        counts = {k: v for k, v in counts.items() if v > 0}
        return max(counts, key=counts.get) if counts else ''

    # ------------------------------------------------------------------
    # Calling
    # ------------------------------------------------------------------

    def call_all(
        self,
        trna_name: str,
        signatures_df: pd.DataFrame,
        pscm_df: Optional[pd.DataFrame] = None,
        reference_seq: Optional[str] = None,
        discover_novel: bool = False,
        min_coverage: int = 50,
        known_mods_df: Optional[pd.DataFrame] = None,
        finalize: bool = True,
    ) -> pd.DataFrame:
        """
        Call modification sites for one tRNA: one row per site.

        Sites with a rate at or above any channel threshold are candidates.
        Each gets a ranked list of compatible modification profiles
        (``candidates``). The top one becomes ``modification`` only when a
        typical position or a specific substitution pattern supports it
        (``identity_support``); otherwise the site is ``novel_candidate``.
        Detected sites are always reported -- identity never removes a
        detection. When *known_mods_df* (MODOMICS, ``linear_position``) is
        given, novel candidates at a known position are relabelled with that
        modification (``source='known_modomics'``).

        With ``finalize=True`` the FDR family is this tRNA. For sample-wide
        FDR (the pipeline), pass ``finalize=False`` for every tRNA, concat,
        and call :meth:`finalize_calls` with the summed :meth:`count_tests`.

        Args:
            trna_name: Name of the tRNA.
            signatures_df: RT signature DataFrame (``position``,
                ``coverage``, ``mismatch_rate`` [substitutions],
                ``gap_rate``, ``rt_stop_pct``).
            pscm_df: Position-Specific Count Matrix DataFrame.
            reference_seq: Reference sequence string.
            discover_novel: Deprecated, ignored. Sites without a supported
                identity are always reported as ``novel_candidate``;
                dropping them would discard detected signal.
            min_coverage: Minimum coverage for a site to be tested.
            known_mods_df: DataFrame with ``linear_position`` and
                ``modification_short_name`` columns from MODOMICS.
            finalize: Apply FDR, channel firing and confidence now.

        Returns:
            DataFrame with ``source`` column ('known', 'known_modomics', or
            'novel_candidate').
        """
        tested = self._test_channels(signatures_df, min_coverage)
        if tested.empty:
            return pd.DataFrame()
        gate_cols = [f'gate_{c}' for c in CHANNELS]
        sites = tested[tested[gate_cols].any(axis=1)]
        if sites.empty:
            return pd.DataFrame()

        correct_nt = (signatures_df.set_index('position')['correct_nt']
                      if 'correct_nt' in signatures_df.columns else None)

        rows = []
        for rec in sites.to_dict('records'):
            position = int(rec['position'])
            ref_nt = None
            if reference_seq is not None and position <= len(reference_seq):
                ref_nt = reference_seq[position - 1].upper()
            elif correct_nt is not None and position in correct_nt.index:
                ref_nt = str(correct_nt.loc[position]).upper()
            pscm_row = None
            if pscm_df is not None and 0 <= position - 1 < len(pscm_df):
                pscm_row = pscm_df.iloc[position - 1]

            gated = {c: bool(rec[f'gate_{c}']) for c in CHANNELS}
            ranked = self._rank_candidates(position, ref_nt, pscm_row, gated,
                                           rec['rate_rt_stop'])
            fraction, typical, support = 0.0, False, 'none'
            if ranked:
                top, fraction, typical = ranked[0]
                specific = (fraction > 0 and bool(top.mismatch_pattern)
                            and not top.mismatch_pattern.endswith('->any'))
                support = ('+'.join(lbl for lbl, ok in
                                    (('position', typical), ('pattern', specific)) if ok)
                           or 'ref_nt_only')
            # A label needs evidence beyond reference-nucleotide compatibility;
            # otherwise the site is detected but its identity is unresolved
            if support in ('none', 'ref_nt_only'):
                modification, full_name, source = (
                    'novel_candidate', 'unknown modification', 'novel_candidate')
            else:
                modification, full_name, source = top.name, top.full_name, 'known'

            row = dict(rec)
            row.update({
                'trna_name': trna_name,
                'modification': modification,
                'full_name': full_name,
                'candidates': ';'.join(p.name for p, _, _ in ranked),
                'n_candidates': len(ranked),
                'identity_support': support,
                'pattern_fraction': fraction,
                'in_typical_position': typical,
                'mismatch_rate': rec['rate_mismatch'],
                'gap_rate': rec['rate_deletion'],
                'rt_stop_pct': rec['rate_rt_stop'] * 100,
                'dominant_pattern': (self._dominant_substitution(ref_nt, pscm_row)
                                     if gated['mismatch'] else ''),
                'source': source,
            })
            rows.append(row)

        if not rows:
            return pd.DataFrame()
        calls = pd.DataFrame(rows)

        # --- MODOMICS-guided relabelling of novel candidates ---
        if (known_mods_df is not None
                and not known_mods_df.empty
                and 'linear_position' in known_mods_df.columns):
            # Build a lookup: linear_position → (short_name, full_name)
            pos_to_mod: Dict[int, Tuple[str, str]] = {}
            for _, row in known_mods_df.iterrows():
                lp = int(row['linear_position'])
                short = row.get('modification_short_name', 'known')
                full = row.get('modification_full_name', short)
                # If multiple mods at same position, concatenate
                if lp in pos_to_mod:
                    prev_short, prev_full = pos_to_mod[lp]
                    if short not in prev_short:
                        pos_to_mod[lp] = (
                            f"{prev_short}; {short}",
                            f"{prev_full}; {full}",
                        )
                else:
                    pos_to_mod[lp] = (short, full)

            novel_mask = calls['source'] == 'novel_candidate'
            for idx in calls.index[novel_mask]:
                pos = int(calls.at[idx, 'position'])
                if pos in pos_to_mod:
                    short_name, full_name = pos_to_mod[pos]
                    calls.at[idx, 'modification'] = short_name
                    calls.at[idx, 'full_name'] = full_name
                    calls.at[idx, 'source'] = 'known_modomics'
                    calls.at[idx, 'in_typical_position'] = True

        if finalize:
            return self.finalize_calls(
                calls, n_tests=self.count_tests(signatures_df, min_coverage))
        return calls

    def finalize_calls(self, calls: pd.DataFrame, n_tests: int) -> pd.DataFrame:
        """Apply per-channel FDR, decide which channels fired, score confidence.

        Args:
            calls: Unfinalised output of :meth:`call_all` (``finalize=False``),
                possibly concatenated across tRNAs of one sample.
            n_tests: Size of the FDR family (summed :meth:`count_tests`).

        Returns:
            Called sites only (at least one channel fired), sorted by
            confidence. Per-channel columns: ``pvalue_*``, ``qvalue_*``,
            ``fold_change_*``; plus ``channels_fired``, ``dominant_channel``,
            ``pvalue`` (Bonferroni over channels, for replicate Fisher
            combination), ``fold_change`` / ``background_error_rate`` of the
            dominant channel, ``confidence``, ``channel_prior_factor`` and
            ``unexpected_channel``.
        """
        if calls.empty:
            return calls
        df = calls.reset_index(drop=True).copy()
        gates = np.column_stack([df[f'gate_{c}'].to_numpy(bool) for c in CHANNELS])

        if self.statistical_test:
            pvals = np.column_stack([df[f'pvalue_{c}'].to_numpy(np.float64)
                                     for c in CHANNELS])
            q = _bh_qvalues(pvals.ravel(), n_tests).reshape(pvals.shape)
            for i, c in enumerate(CHANNELS):
                df[f'qvalue_{c}'] = q[:, i]
            fired = gates & (q < self.alpha)
            df['pvalue'] = np.minimum(1.0, pvals.min(axis=1) * len(CHANNELS))
        else:
            for c in CHANNELS:
                df[f'qvalue_{c}'] = np.nan
            fired = gates
            df['pvalue'] = None

        keep = fired.any(axis=1)
        df, fired = df[keep].reset_index(drop=True), fired[keep]
        if df.empty:
            return pd.DataFrame()

        rates = np.column_stack([df[f'rate_{c}'].to_numpy(np.float64) for c in CHANNELS])
        thresholds = np.array([self.channel_thresholds[c] for c in CHANNELS])
        means = np.array([self.channel_backgrounds[c].mean for c in CHANNELS])

        for i, c in enumerate(CHANNELS):
            df[f'fold_change_{c}'] = (rates[:, i] / means[i] if means[i] > 0 else np.nan)

        # Dominant channel: strongest fired effect relative to its threshold
        rel = np.where(fired, rates / thresholds, -np.inf)
        dom = rel.argmax(axis=1)
        df['channels_fired'] = ['+'.join(c for c, f in zip(CHANNELS, row) if f)
                                for row in fired]
        df['n_channels_fired'] = fired.sum(axis=1)
        df['dominant_channel'] = [CHANNELS[i] for i in dom]
        df['fold_change'] = [rates[j, i] / means[i] if means[i] > 0 else np.nan
                             for j, i in enumerate(dom)]
        df['background_error_rate'] = means[dom]
        df['fdr_significant'] = bool(self.statistical_test)

        # Confidence: ranking score only. Effect size of the strongest fired
        # channel, agreement across channels, position prior, coverage.
        effect = np.where(fired, np.clip((rates - thresholds) / (0.5 - thresholds), 0, 1), 0)
        base = 0.2 + 0.4 * effect.max(axis=1)
        base += 0.1 * np.minimum(fired.sum(axis=1) - 1, 2)
        base += 0.2 * df['in_typical_position'].astype(bool).to_numpy()
        coverage_score = np.minimum(1.0, np.log10(df['coverage'].to_numpy(np.float64) + 1) / 4.0)
        base *= 0.5 + 0.5 * coverage_score

        # Channel priors scale confidence (never detection): signal in the
        # enzyme's usual channel up-weights, a rarely used channel down-weights
        if self.channel_priors:
            weights = np.array([self.channel_priors['weights'][c] for c in CHANNELS])
            fired_w = np.where(fired, weights, 0.0).max(axis=1)
            factor = 1.0 + (fired_w - 1.0 / len(CHANNELS))
            df['unexpected_channel'] = fired_w < UNEXPECTED_CHANNEL_WEIGHT
        else:
            factor = np.ones(len(df))
            df['unexpected_channel'] = False
        df['channel_prior_factor'] = factor
        df['confidence'] = np.minimum(1.0, base * factor)

        df = df[base >= self.min_confidence]
        df = df.drop(columns=[f'{p}_{c}' for p in ('gate', 'rate') for c in CHANNELS])

        # Provenance: which RT produced the data (None if unrecorded)
        df['rt_enzyme'] = self.rt_enzyme
        df['rt_temp'] = self.rt_temp
        return df.sort_values('confidence', ascending=False).reset_index(drop=True)

    def filter_by_confidence(
        self,
        calls_df: pd.DataFrame,
        min_confidence: float
    ) -> pd.DataFrame:
        """
        Filter modification calls by confidence threshold.

        Args:
            calls_df: DataFrame with modification calls
            min_confidence: Minimum confidence threshold

        Returns:
            Filtered DataFrame
        """
        return calls_df[calls_df['confidence'] >= min_confidence].copy()

    def summarize_modifications(
        self,
        calls_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Summarize modification calls by type.

        Args:
            calls_df: DataFrame with modification calls

        Returns:
            Summary DataFrame with counts per modification type
        """
        if calls_df.empty:
            return pd.DataFrame()

        summary = calls_df.groupby('modification').agg({
            'position': 'count',
            'confidence': ['mean', 'std', 'min', 'max'],
            'mismatch_rate': 'mean',
            'coverage': 'mean'
        }).round(3)

        summary.columns = ['_'.join(col).strip('_') for col in summary.columns]
        summary = summary.rename(columns={'position_count': 'num_sites'})
        summary = summary.reset_index()

        return summary


class ReplicateAggregator:
    """Aggregate per-sample modification calls across biological replicates.

    Uses Fisher's combined probability test to merge p-values from
    independent replicate samples and a "double-sieve" filter:
    1. The modification must be detected in >= *min_replicates* samples.
    2. The Fisher combined p-value must be < *alpha*.

    Produces two outputs:
    - **aggregated_modifications**: all sites detected in >=1 replicate,
      with replicate count and Fisher p-value.
    - **consensus_modifications**: the subset passing the double-sieve.
    """

    def __init__(self, min_replicates: int = 3, alpha: float = 0.01):
        self.min_replicates = min_replicates
        self.alpha = alpha

    def aggregate(
        self,
        per_sample_calls: Dict[str, pd.DataFrame],
        replicate_groups: Dict[str, List[str]],
    ) -> pd.DataFrame:
        """Aggregate modification calls across replicate groups.

        Args:
            per_sample_calls: {sample_name_unique: calls_df} — each
                DataFrame has columns including trna_name, position,
                modification, mismatch_rate, fold_change, coverage,
                confidence, pvalue.
            replicate_groups: {condition_name: [sample_name_unique, ...]}.

        Returns:
            DataFrame with aggregated calls (one row per condition x
            trna x position x modification).  Includes
            ``consensus_call`` boolean column.
        """
        rows: List[dict] = []

        for condition, members in replicate_groups.items():
            n_total = len(members)
            # Collect calls from all replicates in this group
            group_dfs = []
            for snu in members:
                df = per_sample_calls.get(snu)
                if df is not None and not df.empty:
                    group_dfs.append(df)

            if not group_dfs:
                continue

            combined = pd.concat(group_dfs, ignore_index=True)

            # Group by modification site
            for (trna, pos, mod), grp in combined.groupby(
                ['trna_name', 'position', 'modification']
            ):
                n_detected = int(grp.shape[0])

                # Fisher combined p-value
                pvals = grp['pvalue'].dropna().values.astype(float)
                pvals = np.clip(pvals, 1e-300, 1.0)
                if len(pvals) >= 2:
                    _, fisher_p = combine_pvalues(pvals, method='fisher')
                elif len(pvals) == 1:
                    fisher_p = float(pvals[0])
                else:
                    fisher_p = np.nan

                mean_mm = float(grp['mismatch_rate'].mean())
                mean_fc = float(grp['fold_change'].mean()) if 'fold_change' in grp.columns else np.nan
                mean_cov = float(grp['coverage'].mean())
                mean_conf = float(grp['confidence'].mean())

                consensus = (
                    n_detected >= self.min_replicates
                    and not np.isnan(fisher_p)
                    and fisher_p < self.alpha
                )

                rows.append({
                    'sample_name': condition,
                    'trna_name': trna,
                    'position': pos,
                    'modification': mod,
                    'n_replicates_detected': n_detected,
                    'n_replicates_total': n_total,
                    'fisher_combined_pvalue': fisher_p,
                    'fisher_significant': (
                        not np.isnan(fisher_p) and fisher_p < self.alpha
                    ),
                    'mean_mismatch_rate': mean_mm,
                    'mean_fold_change': mean_fc,
                    'mean_coverage': mean_cov,
                    'mean_confidence': mean_conf,
                    'consensus_call': consensus,
                })

        if not rows:
            return pd.DataFrame()
        return pd.DataFrame(rows)


def benjamini_hochberg_fdr(pvalues, alpha=0.05):
    """Benjamini-Hochberg FDR correction.

    Manual implementation to avoid scipy version dependency issues.

    Args:
        pvalues: Array-like of p-values.
        alpha: FDR threshold (default 0.05).

    Returns:
        Boolean array indicating which tests pass FDR correction.
    """
    pvals = np.asarray(pvalues, dtype=np.float64)
    n = len(pvals)
    if n == 0:
        return np.array([], dtype=bool)

    # Sort p-values and track original indices
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]

    # BH threshold: p(i) <= (i / n) * alpha
    thresholds = np.arange(1, n + 1) / n * alpha

    # Find largest k where p(k) <= threshold(k)
    below = sorted_pvals <= thresholds
    significant = np.zeros(n, dtype=bool)

    if below.any():
        max_k = np.max(np.where(below)[0])
        # All tests up to and including max_k are significant
        significant[sorted_idx[:max_k + 1]] = True

    return significant
