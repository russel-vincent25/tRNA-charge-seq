"""
Per-RT-enzyme detection-channel priors.

The channel a modification surfaces in (mismatch, deletion, RT stop) is a
property of the reverse transcriptase, not of the modification: on the
2024-06-19 four-RT comparison the same site reads as a deletion on Maxima /
SSIV, a mismatch on TGIRT and an RT stop on Indura.

Priors are *derived* from per-sample ``mismatch_profile`` and ``rt_profile``
tables (raw per-position rates, upstream of any calling logic), never
hand-set. For each enzyme (and enzyme x temperature) they record, among
site-observations carrying signal in any channel, the fraction carrying
signal in each channel; ``weight_*`` normalises those fractions to sum to 1.

Regenerate the shipped table with::

    python -m trnaseq.modifications.channel_priors \\
        /path/to/2024-06-19-RT_comp/modification_analysis

How the weights enter calling is deliberately not decided here; see the
``combined`` aggregation design (CHANNEL_FIX_SPEC change 4).
"""

import argparse
import re
import warnings
from datetime import date
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd

CHANNELS = ('mismatch', 'deletion', 'rt_stop')

# Per-channel signal thresholds, used both for prior derivation and as the
# caller's default calling thresholds. 0.10 in every channel follows the
# published standard: Nakano et al. 2025 (Nat Commun 16:1047) call a
# modification at ">10% of RT misincorporation or stop". On the four-RT data
# an RT-stop threshold of 0.10 rather than 0.20 recovers 14 further sites --
# D16, m7G46, s2C33, Psi40, all RT-stop-driven -- at identical precision by
# both MODOMICS agreement and replicate reproducibility.
DERIVATION_THRESHOLDS = {'mismatch': 0.10, 'deletion': 0.10, 'rt_stop': 0.10}
DERIVATION_MIN_COVERAGE = 100

PRIORS_CSV = Path(__file__).parent / 'data' / 'rt_channel_priors.csv'

# Squashed (lowercase, no space/hyphen/underscore) name -> canonical enzyme
_ENZYME_ALIASES = {
    'maxima': 'Maxima', 'maximah': 'Maxima', 'maximahminus': 'Maxima',
    'ssiv': 'SSIV', 'ss4': 'SSIV', 'superscriptiv': 'SSIV',
    'superscript4': 'SSIV',
    'tgirt': 'TGIRT', 'tgirtiii': 'TGIRT', 'tgirt3': 'TGIRT',
    'indura': 'Indura', 'induro': 'Indura',
}

_COLUMNS = (
    ['rt_enzyme', 'rt_temp', 'n_signal_obs']
    + [f'frac_{c}' for c in CHANNELS]
    + [f'weight_{c}' for c in CHANNELS]
)


def normalize_rt_enzyme(name: Optional[str]) -> Optional[str]:
    """Map an RT enzyme name to its canonical key, or None if unrecognised."""
    if name is None:
        return None
    squashed = re.sub(r'[\s\-_]', '', str(name)).lower()
    return _ENZYME_ALIASES.get(squashed)


def load_profile_dir(
    directory: Union[str, Path],
    sample_pattern: str = r'^(?P<rt_enzyme>[^-]+)-(?P<rt_temp>\d+)-',
    min_coverage: int = DERIVATION_MIN_COVERAGE,
) -> pd.DataFrame:
    """Load per-sample mismatch/RT profiles into one observation table.

    Each sub-directory of *directory* is a sample holding
    ``mismatch_profile.parquet`` and ``rt_profile.parquet``. Enzyme and
    temperature are parsed from the sample name with *sample_pattern*
    (named groups ``rt_enzyme`` and ``rt_temp``).

    Position 1 is dropped: ``rt_stop_fraction`` is 1.0 there by construction.

    Returns:
        One row per (sample, tRNA, position) with ``mismatch``,
        ``deletion`` and ``rt_stop`` rate columns.
    """
    pattern = re.compile(sample_pattern)
    frames = []
    for sample_dir in sorted(Path(directory).iterdir()):
        mm_path = sample_dir / 'mismatch_profile.parquet'
        rt_path = sample_dir / 'rt_profile.parquet'
        if not (mm_path.exists() and rt_path.exists()):
            continue
        m = pattern.match(sample_dir.name)
        if m is None:
            continue
        mm = pd.read_parquet(mm_path)
        rt = pd.read_parquet(rt_path)
        obs = mm.merge(
            rt[['tRNA_name', 'position', 'rt_stop_fraction']],
            on=['tRNA_name', 'position'],
        )
        obs = obs[(obs['coverage'] >= min_coverage) & (obs['position'] > 1)]
        frames.append(pd.DataFrame({
            'sample': sample_dir.name,
            'rt_enzyme': normalize_rt_enzyme(m.group('rt_enzyme')),
            'rt_temp': float(m.group('rt_temp')),
            'tRNA_name': obs['tRNA_name'].values,
            'position': obs['position'].values,
            'coverage': obs['coverage'].values,
            'mismatch': obs['mismatch_rate'].values,
            'deletion': obs['deletion_rate'].values,
            'rt_stop': obs['rt_stop_fraction'].values,
        }))
    if not frames:
        raise FileNotFoundError(f"No sample profiles found under {directory}")
    return pd.concat(frames, ignore_index=True)


def derive_channel_priors(
    obs_df: pd.DataFrame,
    thresholds: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """Derive per-enzyme channel priors from an observation table.

    Args:
        obs_df: Output of :func:`load_profile_dir` (columns ``rt_enzyme``,
            ``rt_temp`` and one rate column per channel).
        thresholds: Per-channel signal thresholds; defaults to
            :data:`DERIVATION_THRESHOLDS`.

    Returns:
        One row per enzyme x temperature plus one pooled row per enzyme
        (``rt_temp`` NaN).
    """
    thresholds = thresholds or DERIVATION_THRESHOLDS
    hits = pd.DataFrame({c: obs_df[c] >= thresholds[c] for c in CHANNELS})
    hits['rt_enzyme'] = obs_df['rt_enzyme']
    hits['rt_temp'] = obs_df['rt_temp']
    signal = hits[hits[list(CHANNELS)].any(axis=1)]

    rows = []
    groups = [(k, g) for k, g in signal.groupby(['rt_enzyme', 'rt_temp'])]
    groups += [((k, np.nan), g) for k, g in signal.groupby('rt_enzyme')]
    for (enzyme, temp), g in groups:
        fracs = {c: float(g[c].mean()) for c in CHANNELS}
        total = sum(fracs.values())
        row = {'rt_enzyme': enzyme, 'rt_temp': temp, 'n_signal_obs': len(g)}
        row.update({f'frac_{c}': fracs[c] for c in CHANNELS})
        row.update({f'weight_{c}': fracs[c] / total for c in CHANNELS})
        rows.append(row)

    priors = pd.DataFrame(rows, columns=_COLUMNS)
    return priors.sort_values(
        ['rt_enzyme', 'rt_temp'], na_position='first'
    ).reset_index(drop=True)


def write_channel_priors(
    priors: pd.DataFrame,
    path: Union[str, Path] = PRIORS_CSV,
    source: str = '',
) -> None:
    """Write priors as CSV with a provenance comment header."""
    header = [
        '# Per-RT channel priors derived by trnaseq.modifications.channel_priors',
        f'# generated: {date.today().isoformat()}',
        f'# source: {source}',
        '# thresholds: ' + ', '.join(
            f'{c}>={t}' for c, t in DERIVATION_THRESHOLDS.items()),
        f'# min_coverage: {DERIVATION_MIN_COVERAGE}; position 1 excluded',
        '# frac_* = fraction of signal site-observations with signal in that '
        'channel (rows may exceed 1); weight_* = frac_* normalised to sum 1',
        '# rt_temp empty = pooled over temperatures',
    ]
    with open(path, 'w') as fh:
        fh.write('\n'.join(header) + '\n')
        priors.round(4).to_csv(fh, index=False)


@lru_cache(maxsize=4)
def _read_priors(path: str) -> pd.DataFrame:
    return pd.read_csv(path, comment='#')


def load_channel_priors(
    rt_enzyme: Optional[str],
    rt_temp: Optional[float] = None,
    path: Union[str, Path] = PRIORS_CSV,
) -> Optional[Dict]:
    """Look up channel priors for an enzyme (and temperature).

    Uses the enzyme x temperature row when *rt_temp* matches a derived
    temperature, otherwise the enzyme's pooled row.

    Returns:
        ``{'rt_enzyme', 'rt_temp', 'n_signal_obs', 'weights': {channel: w}}``,
        or ``None`` when the enzyme is unset or has no derived priors --
        callers then fall back to enzyme-agnostic 'combined' behaviour.
    """
    if rt_enzyme is None:
        return None
    canonical = normalize_rt_enzyme(rt_enzyme)
    priors = _read_priors(str(path))
    rows = priors[priors['rt_enzyme'] == canonical] if canonical else priors.iloc[:0]
    if rows.empty:
        warnings.warn(
            f"No channel priors for RT enzyme '{rt_enzyme}'; known: "
            f"{sorted(priors['rt_enzyme'].unique())}. Falling back to "
            f"enzyme-agnostic 'combined' calling."
        )
        return None

    row = None
    if rt_temp is not None:
        match = rows[np.isclose(rows['rt_temp'], float(rt_temp))]
        if not match.empty:
            row = match.iloc[0]
    if row is None:
        row = rows[rows['rt_temp'].isna()].iloc[0]

    return {
        'rt_enzyme': canonical,
        'rt_temp': None if pd.isna(row['rt_temp']) else float(row['rt_temp']),
        'n_signal_obs': int(row['n_signal_obs']),
        'weights': {c: float(row[f'weight_{c}']) for c in CHANNELS},
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Derive per-RT channel priors from per-sample '
                    'mismatch_profile/rt_profile parquets.')
    parser.add_argument('profile_dir',
                        help='Directory of <enzyme>-<temp>-<rep>/ sample dirs')
    parser.add_argument('-o', '--output', default=str(PRIORS_CSV))
    args = parser.parse_args(argv)

    obs = load_profile_dir(args.profile_dir)
    priors = derive_channel_priors(obs)
    write_channel_priors(priors, args.output, source=args.profile_dir)
    print(f"{len(obs):,} site-observations -> {args.output}")
    print(priors.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
