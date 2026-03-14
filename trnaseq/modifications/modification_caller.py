"""
Modification Caller for tRNA RT Signatures

This module annotates RT signatures with specific tRNA modification types based on
known RT signature patterns from the literature.

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
from dataclasses import dataclass
from scipy.stats import binomtest, combine_pvalues


@dataclass
class ModificationProfile:
    """
    Profile for a specific tRNA modification type.

    Attributes:
        name: Modification name (e.g., 'm1A')
        full_name: Full chemical name
        typical_positions: Common positions where this modification occurs
        signature_type: Type of RT signature ('mismatch', 'rt_stop', 'gap', 'combined')
        mismatch_pattern: Expected mismatch pattern (e.g., 'A->G', 'A->any')
        min_rate: Minimum rate to call this modification
        rt_stop_required: Whether RT stops are required for calling
        min_rt_stop_pct: Minimum RT stop percentage if required
    """
    name: str
    full_name: str
    typical_positions: List[int]
    signature_type: str
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
        signature_type='combined',
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
        signature_type='mismatch',
        mismatch_pattern='C->T',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.2
    ),

    'pseudouridine': ModificationProfile(
        name='Ψ',
        full_name='pseudouridine',
        typical_positions=[27, 28, 31, 32, 39, 40, 55, 13, 38],
        signature_type='mismatch',
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm7G': ModificationProfile(
        name='m7G',
        full_name='7-methylguanosine',
        typical_positions=[46],
        signature_type='mismatch',
        mismatch_pattern='G->any',  # TGIRT/Maxima produce G->C; some enzymes G->A
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.1
    ),

    'm5C': ModificationProfile(
        name='m5C',
        full_name='5-methylcytosine',
        typical_positions=[48, 49, 34, 40],
        signature_type='mismatch',
        mismatch_pattern='C->T',
        min_rate=0.05,  # Subtle signature
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'i6A': ModificationProfile(
        name='i6A',
        full_name='N6-isopentenyladenosine',
        typical_positions=[37],
        signature_type='mismatch',
        mismatch_pattern='A->G',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm2G': ModificationProfile(
        name='m2G',
        full_name='N2-methylguanosine',
        typical_positions=[10, 26],
        signature_type='mismatch',
        mismatch_pattern='G->A',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    'm22G': ModificationProfile(
        name='m22G',
        full_name='N2,N2-dimethylguanosine',
        typical_positions=[26],
        signature_type='mismatch',
        mismatch_pattern='G->A',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    's4U': ModificationProfile(
        name='s4U',
        full_name='4-thiouridine',
        typical_positions=[8, 9, 4],
        signature_type='mismatch',
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm1G': ModificationProfile(
        name='m1G',
        full_name='1-methylguanosine',
        typical_positions=[37, 9],
        signature_type='combined',
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
        signature_type='mismatch',
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    'mnm5s2U': ModificationProfile(
        name='mnm5s2U',
        full_name='5-methylaminomethyl-2-thiouridine',
        typical_positions=[34],
        signature_type='combined',
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    't6A': ModificationProfile(
        name='t6A',
        full_name='N6-threonylcarbamoyladenosine',
        typical_positions=[37],
        signature_type='mismatch',
        mismatch_pattern='A->T',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'ms2i6A': ModificationProfile(
        name='ms2i6A',
        full_name='2-methylthio-N6-isopentenyladenosine',
        typical_positions=[37],
        signature_type='combined',
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
        signature_type='mismatch',
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
        signature_type='mismatch',
        mismatch_pattern='C->T',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    'Gm': ModificationProfile(
        name='Gm',
        full_name="2'-O-methylguanosine",
        typical_positions=[18, 34],
        signature_type='mismatch',
        mismatch_pattern='G->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'Cm': ModificationProfile(
        name='Cm',
        full_name="2'-O-methylcytidine",
        typical_positions=[32, 34],
        signature_type='mismatch',
        mismatch_pattern='C->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'Um': ModificationProfile(
        name='Um',
        full_name="2'-O-methyluridine",
        typical_positions=[32, 44],
        signature_type='mismatch',
        mismatch_pattern='U->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'Am': ModificationProfile(
        name='Am',
        full_name="2'-O-methyladenosine",
        typical_positions=[44],
        signature_type='mismatch',
        mismatch_pattern='A->any',
        min_rate=0.05,
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'Q': ModificationProfile(
        name='Q',
        full_name='queuosine',
        typical_positions=[34],
        signature_type='mismatch',
        mismatch_pattern='G->any',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm6A': ModificationProfile(
        name='m6A',
        full_name='N6-methyladenosine',
        typical_positions=[37, 58],
        signature_type='mismatch',
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


class ModificationCaller:
    """
    Call specific tRNA modifications based on RT signature patterns.

    This class takes RT signature data and matches it against known modification
    profiles to identify likely modification sites with confidence scores.

    Example:
        >>> from trnaseq.modifications import RTSignatureAnalyzer, ModificationCaller
        >>>
        >>> # Get RT signatures
        >>> analyzer = RTSignatureAnalyzer()
        >>> results = analyzer.analyze_all_trnas()
        >>>
        >>> # Call modifications
        >>> caller = ModificationCaller(organism='human')
        >>> modifications = caller.call_modifications_for_trna(
        ...     'tRNA-Thr-AGT-1-1',
        ...     results['tRNA-Thr-AGT-1-1']['signatures']
        ... )
    """

    def __init__(
        self,
        organism: str = 'human',
        min_confidence: float = 0.5,
        use_position_priors: bool = True,
        statistical_test: bool = True,
        alpha: float = 0.01,
        background_error_rate: float = 0.01,
    ):
        """
        Initialize modification caller.

        Args:
            organism: Organism name (for position priors)
            min_confidence: Minimum confidence score to report (0-1)
            use_position_priors: Use known modification positions to boost confidence
            statistical_test: Perform binomial test for significance
            alpha: Significance level for statistical test
            background_error_rate: Expected sequencing error rate for binomial
                test and fold-change computation. Use
                :func:`estimate_background_error_rate` to set empirically.
        """
        self.organism = organism
        self.min_confidence = min_confidence
        self.use_position_priors = use_position_priors
        self.statistical_test = statistical_test
        self.alpha = alpha
        self.background_error_rate = background_error_rate

        # Load modification profiles
        self.profiles = MODIFICATION_PROFILES

    def match_mismatch_pattern(
        self,
        pattern: str,
        ref_nt: str,
        pscm_row: pd.Series
    ) -> Tuple[bool, float]:
        """
        Check if observed mismatch pattern matches expected pattern.

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
        if expected_ref != ref_nt and not (
            {expected_ref, ref_nt} <= {'U', 'T'}
        ):
            return False, 0.0

        # Calculate total mismatches
        total_coverage = pscm_row.sum()
        if total_coverage == 0:
            return False, 0.0

        correct_count = pscm_row.get(ref_nt, 0)
        total_mismatches = total_coverage - correct_count

        if total_mismatches == 0:
            return False, 0.0

        # Check if specific nucleotide pattern matches
        if expected_obs == 'any':
            # Any mismatch is acceptable
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
            # Pattern matches if this specific mismatch is dominant
            matches = obs_count >= (total_mismatches * 0.5)  # At least 50% of mismatches
            return matches, fraction

    def calculate_confidence(
        self,
        profile: ModificationProfile,
        position: int,
        mismatch_rate: float,
        rt_stop_pct: float,
        gap_rate: float,
        pattern_fraction: float,
        coverage: float
    ) -> float:
        """
        Calculate confidence score for modification call.

        Confidence is based on:
        1. Mismatch rate (higher is better)
        2. RT stop percentage (if required)
        3. Pattern specificity
        4. Position match (if using priors)
        5. Coverage (higher is more reliable)

        Args:
            profile: Modification profile
            position: Position in tRNA
            mismatch_rate: Observed mismatch rate
            rt_stop_pct: Observed RT stop percentage
            gap_rate: Observed gap rate
            pattern_fraction: Fraction matching expected pattern
            coverage: Read coverage

        Returns:
            Confidence score (0-1)
        """
        confidence = 0.0

        # Base confidence from mismatch rate
        if mismatch_rate >= profile.min_rate:
            # Scale between min_rate and 0.5 (saturate at 50% mismatch)
            confidence += min(1.0, (mismatch_rate - profile.min_rate) / (0.5 - profile.min_rate)) * 0.4

        # RT stop contribution
        if profile.rt_stop_required:
            if rt_stop_pct >= profile.min_rt_stop_pct:
                rt_stop_score = min(1.0, rt_stop_pct / 50.0)  # Saturate at 50%
                confidence += rt_stop_score * 0.3
            else:
                # Penalize if RT stop required but not present
                confidence *= 0.3
        elif profile.signature_type == 'combined' and rt_stop_pct >= profile.min_rt_stop_pct:
            # Bonus for combined-type profiles when RT stop is present
            rt_stop_score = min(1.0, rt_stop_pct / 50.0)
            confidence += rt_stop_score * 0.2

        # Pattern specificity
        if pattern_fraction > 0:
            confidence += pattern_fraction * 0.2

        # Position prior
        if self.use_position_priors and position in profile.typical_positions:
            confidence += 0.2

        # Coverage contribution (higher coverage = more reliable)
        coverage_score = min(1.0, np.log10(coverage + 1) / 4.0)  # Saturate at 10,000 reads
        confidence *= (0.5 + 0.5 * coverage_score)  # Scale between 50-100% based on coverage

        # Apply profile-specific weight
        confidence *= profile.confidence_weight

        # Cap at 1.0
        return min(1.0, confidence)

    def perform_statistical_test(
        self,
        coverage: int,
        mismatch_count: int,
        expected_error_rate: float = None,
    ) -> float:
        """
        Perform binomial test to check if mismatch rate is significantly elevated.

        Args:
            coverage: Total read coverage
            mismatch_count: Number of mismatches observed
            expected_error_rate: Expected sequencing error rate. If *None*,
                uses ``self.background_error_rate``.

        Returns:
            P-value from binomial test
        """
        if expected_error_rate is None:
            expected_error_rate = self.background_error_rate
        if coverage == 0:
            return 1.0

        # Binomial test: is mismatch rate significantly > error rate?
        result = binomtest(
            k=int(mismatch_count),
            n=int(coverage),
            p=expected_error_rate,
            alternative='greater'
        )

        return result.pvalue

    def call_modification_at_position(
        self,
        position: int,
        signature_row: pd.Series,
        pscm_row: Optional[pd.Series] = None,
        ref_nt: Optional[str] = None
    ) -> List[Dict]:
        """
        Call modifications at a specific position.

        Args:
            position: Position in tRNA (1-based)
            signature_row: Row from signature DataFrame
            pscm_row: Row from PSCM (for pattern matching)
            ref_nt: Reference nucleotide at this position

        Returns:
            List of modification calls with confidence scores
        """
        calls = []

        mismatch_rate = signature_row.get('mismatch_rate', 0)
        rt_stop_pct = signature_row.get('rt_stop_pct', 0)
        gap_rate = signature_row.get('gap_rate', 0)
        coverage = signature_row.get('coverage', 0)

        # Try each modification profile
        for mod_name, profile in self.profiles.items():
            # Check if mismatch rate meets minimum
            if mismatch_rate < profile.min_rate:
                continue

            # Check if RT stop is required
            if profile.rt_stop_required and rt_stop_pct < profile.min_rt_stop_pct:
                continue

            # Check mismatch pattern if PSCM provided
            pattern_matches = True
            pattern_fraction = 0.0
            if pscm_row is not None and ref_nt is not None and profile.mismatch_pattern:
                pattern_matches, pattern_fraction = self.match_mismatch_pattern(
                    profile.mismatch_pattern,
                    ref_nt,
                    pscm_row
                )
                if not pattern_matches:
                    continue

            # Calculate confidence
            confidence = self.calculate_confidence(
                profile,
                position,
                mismatch_rate,
                rt_stop_pct,
                gap_rate,
                pattern_fraction,
                coverage
            )

            # Statistical test if requested
            pvalue = None
            if self.statistical_test and coverage > 0:
                mismatch_count = mismatch_rate * coverage
                pvalue = self.perform_statistical_test(coverage, mismatch_count)

                # Reduce confidence if not significant
                if pvalue > self.alpha:
                    confidence *= 0.5

            # Only report if meets minimum confidence
            if confidence >= self.min_confidence:
                call = {
                    'position': position,
                    'modification': profile.name,
                    'full_name': profile.full_name,
                    'confidence': confidence,
                    'mismatch_rate': mismatch_rate,
                    'rt_stop_pct': rt_stop_pct,
                    'gap_rate': gap_rate,
                    'coverage': coverage,
                    'pattern_fraction': pattern_fraction,
                    'pvalue': pvalue,
                    'in_typical_position': position in profile.typical_positions,
                    'fold_change': (mismatch_rate / self.background_error_rate
                                    if self.background_error_rate > 0
                                    else np.nan),
                    'background_error_rate': self.background_error_rate,
                }
                calls.append(call)

        return calls

    def call_modifications_for_trna(
        self,
        trna_name: str,
        signatures_df: pd.DataFrame,
        pscm_df: Optional[pd.DataFrame] = None,
        reference_seq: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Call modifications for a single tRNA.

        Args:
            trna_name: Name of the tRNA
            signatures_df: RT signature DataFrame from RTSignatureAnalyzer
            pscm_df: Position-Specific Count Matrix (optional, for pattern matching)
            reference_seq: Reference sequence (optional, for pattern matching)

        Returns:
            DataFrame with modification calls
        """
        all_calls = []

        # Process each position with a signature
        for _, row in signatures_df.iterrows():
            if not row.get('has_signature', False):
                continue

            position = int(row['position'])

            # Get PSCM row and reference nt if available
            pscm_row = None
            ref_nt = None
            if pscm_df is not None:
                pscm_row = pscm_df.iloc[position - 1]  # Convert to 0-based
            if reference_seq is not None and position <= len(reference_seq):
                ref_nt = reference_seq[position - 1]  # Convert to 0-based

            # Call modifications at this position
            position_calls = self.call_modification_at_position(
                position, row, pscm_row, ref_nt
            )

            for call in position_calls:
                call['trna_name'] = trna_name
                all_calls.append(call)

        if not all_calls:
            return pd.DataFrame()

        # Convert to DataFrame and sort by confidence
        calls_df = pd.DataFrame(all_calls)
        calls_df = calls_df.sort_values(['confidence', 'mismatch_rate'], ascending=False)

        return calls_df

    def call_modifications_for_all_trnas(
        self,
        rt_signature_results: Dict[str, Dict[str, pd.DataFrame]],
        pscm_dict: Optional[Dict[str, pd.DataFrame]] = None,
        reference_sequences: Optional[Dict[str, str]] = None
    ) -> pd.DataFrame:
        """
        Call modifications for all tRNAs.

        Args:
            rt_signature_results: Results from RTSignatureAnalyzer.analyze_all_trnas()
            pscm_dict: Dictionary of tRNA -> PSCM DataFrame
            reference_sequences: Dictionary of tRNA -> sequence

        Returns:
            Combined DataFrame with all modification calls
        """
        all_calls = []

        for trna_name, results in rt_signature_results.items():
            signatures_df = results['signatures']

            pscm_df = pscm_dict.get(trna_name) if pscm_dict else None
            ref_seq = reference_sequences.get(trna_name) if reference_sequences else None

            trna_calls = self.call_modifications_for_trna(
                trna_name,
                signatures_df,
                pscm_df,
                ref_seq
            )

            if not trna_calls.empty:
                all_calls.append(trna_calls)

        if not all_calls:
            return pd.DataFrame()

        combined_df = pd.concat(all_calls, ignore_index=True)
        return combined_df

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

    def call_novel_positions(
        self,
        trna_name: str,
        signatures_df: pd.DataFrame,
        pscm_df: Optional[pd.DataFrame] = None,
        reference_seq: Optional[str] = None,
        min_coverage: int = 50,
    ) -> pd.DataFrame:
        """
        Identify novel modification candidates at signature positions.

        Iterates positions where ``has_signature=True`` and skips those that
        already matched a known profile via :meth:`call_modification_at_position`.
        Remaining positions with statistically significant mismatch rates
        are reported as ``novel_candidate``.

        Args:
            trna_name: Name of the tRNA.
            signatures_df: RT signature DataFrame (must have ``has_signature``).
            pscm_df: Position-Specific Count Matrix DataFrame.
            reference_seq: Reference sequence string.
            min_coverage: Minimum read coverage to consider a position.

        Returns:
            DataFrame of novel modification candidates with dominant mutation
            pattern and statistical evidence.
        """
        all_calls = []

        for _, row in signatures_df.iterrows():
            if not row.get('has_signature', False):
                continue

            position = int(row['position'])
            coverage = row.get('coverage', 0)
            if coverage < min_coverage:
                continue

            # Check if any known profile matches
            pscm_row = None
            ref_nt = None
            if pscm_df is not None:
                pos_idx = position - 1
                if 0 <= pos_idx < len(pscm_df):
                    pscm_row = pscm_df.iloc[pos_idx]
            if reference_seq is not None and position <= len(reference_seq):
                ref_nt = reference_seq[position - 1]

            known_calls = self.call_modification_at_position(
                position, row, pscm_row, ref_nt
            )
            if known_calls:
                continue  # Already matched a known profile

            # Statistical test
            mismatch_rate = row.get('mismatch_rate', 0)
            mismatch_count = mismatch_rate * coverage
            pvalue = self.perform_statistical_test(
                int(coverage), int(mismatch_count)
            )
            if pvalue > self.alpha:
                continue  # Not significant

            # Extract dominant mutation pattern
            dominant_pattern = ''
            if pscm_row is not None and ref_nt is not None:
                nt_counts = {}
                for nt in ['A', 'C', 'G', 'T']:
                    if nt != ref_nt:
                        cnt = pscm_row.get(nt, 0)
                        if cnt > 0:
                            nt_counts[f'{ref_nt}->{nt}'] = cnt
                if nt_counts:
                    dominant_pattern = max(nt_counts, key=nt_counts.get)

            all_calls.append({
                'trna_name': trna_name,
                'position': position,
                'modification': 'novel_candidate',
                'full_name': 'unknown modification',
                'confidence': min(1.0, mismatch_rate * 2),
                'mismatch_rate': mismatch_rate,
                'rt_stop_pct': row.get('rt_stop_pct', 0),
                'gap_rate': row.get('gap_rate', 0),
                'coverage': coverage,
                'pattern_fraction': 0.0,
                'pvalue': pvalue,
                'in_typical_position': False,
                'dominant_pattern': dominant_pattern,
                'source': 'novel_candidate',
                'fold_change': (mismatch_rate / self.background_error_rate
                                if self.background_error_rate > 0
                                else np.nan),
                'background_error_rate': self.background_error_rate,
            })

        return pd.DataFrame(all_calls) if all_calls else pd.DataFrame()

    def call_all(
        self,
        trna_name: str,
        signatures_df: pd.DataFrame,
        pscm_df: Optional[pd.DataFrame] = None,
        reference_seq: Optional[str] = None,
        discover_novel: bool = False,
        min_coverage: int = 50,
        known_mods_df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """
        Run both known profile matching and optional novel discovery.

        Wrapper that calls :meth:`call_modifications_for_trna` for known
        profiles and :meth:`call_novel_positions` for novel candidates,
        then returns the combined results.

        When *known_mods_df* is provided (a DataFrame of MODOMICS-derived
        modifications with a ``linear_position`` column), novel candidates
        whose position matches a known modification are relabelled with the
        known modification name and their confidence is boosted.

        Args:
            trna_name: Name of the tRNA.
            signatures_df: RT signature DataFrame.
            pscm_df: Position-Specific Count Matrix DataFrame.
            reference_seq: Reference sequence string.
            discover_novel: If True, also run novel modification discovery.
            min_coverage: Minimum coverage for novel discovery.
            known_mods_df: DataFrame with ``linear_position`` and
                ``modification_short_name`` columns from MODOMICS.

        Returns:
            Combined DataFrame with ``source`` column ('known',
            'known_modomics', or 'novel_candidate').
        """
        known_df = self.call_modifications_for_trna(
            trna_name, signatures_df, pscm_df, reference_seq
        )
        if not known_df.empty:
            known_df['source'] = 'known'

        if discover_novel:
            novel_df = self.call_novel_positions(
                trna_name, signatures_df, pscm_df, reference_seq,
                min_coverage=min_coverage,
            )
            if not novel_df.empty and not known_df.empty:
                combined = pd.concat([known_df, novel_df], ignore_index=True)
            elif not novel_df.empty:
                combined = novel_df
            else:
                combined = known_df
        else:
            combined = known_df

        # --- MODOMICS-guided relabelling of novel candidates ---
        if (not combined.empty
                and known_mods_df is not None
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

            novel_mask = combined['source'] == 'novel_candidate'
            for idx in combined.index[novel_mask]:
                pos = int(combined.at[idx, 'position'])
                if pos in pos_to_mod:
                    short_name, full_name = pos_to_mod[pos]
                    combined.at[idx, 'modification'] = short_name
                    combined.at[idx, 'full_name'] = full_name
                    combined.at[idx, 'source'] = 'known_modomics'
                    # Boost confidence for MODOMICS-confirmed positions
                    combined.at[idx, 'confidence'] = min(
                        1.0, combined.at[idx, 'confidence'] * 1.5
                    )
                    combined.at[idx, 'in_typical_position'] = True

        if not combined.empty:
            # Apply FDR correction across all tested positions
            if 'pvalue' in combined.columns:
                pvals = combined['pvalue'].fillna(1.0).values
                combined['fdr_significant'] = benjamini_hochberg_fdr(
                    pvals, alpha=self.alpha
                )

        return combined

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
