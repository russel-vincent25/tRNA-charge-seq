"""
Novel-site summary: what the unexplained modification calls look like.

A ``novel_candidate`` call is a detected site with no MODOMICS modification
at the base its signal implicates and no matching substitution pattern.
Per-sample call tables scatter each such site across samples, so this
module collapses them to one row per (tRNA, position) and records what
would help decide what it is:

- which channels fired, and how often (``channels``)
- the substitution spectrum, from ``dominant_pattern`` (``substitutions``)
- median and maximum rate per channel
- which samples, conditions and RT enzymes show it
- whether the same site is labelled in other samples (``labelled_elsewhere``)
- the nearest MODOMICS modification and its offset (``nearest_known``,
  ``nearest_known_offset``), since RT signal often lands a base or two 3'
  of the modified nucleotide

Output of :func:`summarize_novel_sites` is written by stage 6 as
``results/modifications/novel_sites`` and shown in the modification report.
"""

from collections import Counter
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

NOVEL = 'novel_candidate'

_COLUMNS = [
    'trna_name', 'position', 'modified_position', 'ref_nt',
    'n_observations', 'n_samples', 'n_conditions', 'rt_enzymes',
    'dominant_channel', 'channels', 'substitutions',
    'median_mismatch_rate', 'max_mismatch_rate',
    'median_gap_rate', 'max_gap_rate',
    'median_rt_stop_pct', 'max_rt_stop_pct',
    'median_coverage', 'candidates', 'labelled_elsewhere',
    'nearest_known', 'nearest_known_position', 'nearest_known_offset',
    'samples', 'conditions',
]


def _counts(values) -> str:
    """'a:3; b:1' for non-empty values, most frequent first."""
    c = Counter(v for v in values if isinstance(v, str) and v)
    return '; '.join(f'{k}:{n}' for k, n in c.most_common())


def _known_positions(known: Optional[pd.DataFrame]) -> Dict[int, str]:
    if known is None or known.empty or 'linear_position' not in known.columns:
        return {}
    out: Dict[int, List[str]] = {}
    for pos, name in zip(known['linear_position'], known['modification_short_name']):
        out.setdefault(int(pos), [])
        if name not in out[int(pos)]:
            out[int(pos)].append(name)
    return {p: '; '.join(n) for p, n in out.items()}


def summarize_novel_sites(
    per_sample_calls: Dict[str, pd.DataFrame],
    replicate_groups: Optional[Dict[str, List[str]]] = None,
    known_mods: Optional[Dict[str, pd.DataFrame]] = None,
    ref_dict: Optional[Dict[str, dict]] = None,
    window: int = 3,
) -> pd.DataFrame:
    """Collapse ``novel_candidate`` calls to one row per site.

    Args:
        per_sample_calls: {sample_name_unique: calls DataFrame} from
            :meth:`ModificationCaller.finalize_calls`.
        replicate_groups: {condition: [sample_name_unique, ...]}.
        known_mods: {trna_name: MODOMICS map from
            ``MODOMICSAnnotator.get_known_mods_linear``}.
        ref_dict: {trna_name: {'seq': str}} for the reference base.
        window: Report the nearest known modification within this many
            bases of ``modified_position``.

    Returns:
        DataFrame with the columns listed in ``_COLUMNS``, sorted by
        ``n_observations`` (descending). Empty if there are no novel calls.
    """
    frames = []
    for sample, df in per_sample_calls.items():
        if df is None or df.empty or 'source' not in df.columns:
            continue
        frames.append(df.assign(_sample=sample))
    if not frames:
        return pd.DataFrame(columns=_COLUMNS)
    calls = pd.concat(frames, ignore_index=True)

    novel = calls[calls['source'] == NOVEL]
    if novel.empty:
        return pd.DataFrame(columns=_COLUMNS)
    labelled = calls[calls['source'] != NOVEL]
    labelled_at = {k: _counts(g['modification'])
                   for k, g in labelled.groupby(['trna_name', 'position'])}

    condition_of = {s: c for c, members in (replicate_groups or {}).items() for s in members}
    known_mods = known_mods or {}
    ref_dict = ref_dict or {}

    def col(g, name):
        return g[name] if name in g.columns else pd.Series(np.nan, index=g.index)

    rows = []
    for (trna, pos), g in novel.groupby(['trna_name', 'position']):
        pos = int(pos)
        mod_pos = int(col(g, 'modified_position').fillna(pos).mode().iloc[0])
        seq = ref_dict.get(trna, {}).get('seq')
        samples = sorted(set(g['_sample']))
        conditions = sorted({condition_of.get(s, s) for s in samples})

        known = _known_positions(known_mods.get(trna))
        near_name, near_pos, near_off = '', np.nan, np.nan
        if known:
            nearest = min(known, key=lambda p: (abs(p - mod_pos), p))
            if abs(nearest - mod_pos) <= window:
                near_name, near_pos, near_off = known[nearest], nearest, nearest - mod_pos

        rows.append({
            'trna_name': trna,
            'position': pos,
            'modified_position': mod_pos,
            'ref_nt': seq[pos - 1].upper() if seq and pos <= len(seq) else '',
            'n_observations': len(g),
            'n_samples': len(samples),
            'n_conditions': len(conditions),
            'rt_enzymes': ';'.join(sorted({str(e) for e in col(g, 'rt_enzyme').dropna()})),
            'dominant_channel': col(g, 'dominant_channel').mode().iloc[0]
            if col(g, 'dominant_channel').notna().any() else '',
            'channels': _counts(col(g, 'channels_fired')),
            'substitutions': _counts(col(g, 'dominant_pattern')),
            'median_mismatch_rate': g['mismatch_rate'].median(),
            'max_mismatch_rate': g['mismatch_rate'].max(),
            'median_gap_rate': g['gap_rate'].median(),
            'max_gap_rate': g['gap_rate'].max(),
            'median_rt_stop_pct': g['rt_stop_pct'].median(),
            'max_rt_stop_pct': g['rt_stop_pct'].max(),
            'median_coverage': g['coverage'].median(),
            'candidates': _counts(c for cs in col(g, 'candidates').dropna()
                                  for c in str(cs).split(';')),
            'labelled_elsewhere': labelled_at.get((trna, pos), ''),
            'nearest_known': near_name,
            'nearest_known_position': near_pos,
            'nearest_known_offset': near_off,
            'samples': ';'.join(samples),
            'conditions': ';'.join(conditions),
        })

    out = pd.DataFrame(rows, columns=_COLUMNS)
    return out.sort_values(['n_observations', 'trna_name', 'position'],
                           ascending=[False, True, True]).reset_index(drop=True)
