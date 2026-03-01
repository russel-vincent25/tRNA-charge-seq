"""
FragmentAnalyser — per-tRNA fragment classification, RT drop-off profiling,
and fragment length distributions.

Reads per-sample stats CSVs (*_stats.csv.bz2) from the stats_collection
directory and produces:

1. fragment_counts   — per-tRNA per-sample fragment type breakdown
2. rt_dropoff_positions — per-tRNA per-position RT stop histogram
3. fragment_lengths  — fragment length distributions by type
4. fragment_summary  — sample-level summary statistics
"""

from pathlib import Path
import pandas as pd
import numpy as np


# Columns we actually need from the 28-column per-read stats CSV
_USECOLS = [
    'sample_name_unique', 'tRNA_annotation', 'tRNA_annotation_len',
    'amino_acid', 'codon', 'anticodon',
    'align_5p_idx', 'align_3p_idx',
    '5p_cover', '3p_cover',
    'count', 'unique_annotation',
]

# Fragment type classification based on 5p_cover x 3p_cover
_FRAG_TYPES = {
    (True, True): 'full_length',
    (False, True): 'rt_dropoff',
    (True, False): '5p_tRF',
    (False, False): 'degraded',
}


class FragmentAnalyser:
    """Classify tRNA-seq reads into fragment types and profile RT drop-off.

    Parameters
    ----------
    stats_dir : str or Path
        Directory containing ``{sample}_stats.csv.bz2`` files.
    sample_names : list of str, optional
        Samples to process.  If None, auto-detect from filenames.
    count_col : str
        Column to use as read count (default ``'count'``).
    min_reads : int
        Minimum total reads per tRNA to include in output (default 10).
    """

    def __init__(self, stats_dir, sample_names=None, count_col='count',
                 min_reads=10):
        self.stats_dir = Path(stats_dir)
        self.count_col = count_col
        self.min_reads = min_reads

        if sample_names is not None:
            self.sample_names = list(sample_names)
        else:
            self.sample_names = self._discover_samples()

        # Result containers (populated by run())
        self._fragment_counts = None
        self._rt_dropoff = None
        self._fragment_lengths = None
        self._summary = None
        self._coverage = None

    def _discover_samples(self):
        """Auto-detect sample names from *_stats.csv.bz2 files."""
        names = []
        for p in sorted(self.stats_dir.glob('*_stats.csv.bz2')):
            name = p.name.replace('_stats.csv.bz2', '')
            if name != 'ALL_stats_aggregate':
                names.append(name)
        return names

    def run(self):
        """Iterate over samples and populate all result DataFrames."""
        frag_parts = []
        rt_parts = []
        len_parts = []
        cov_parts = []
        summary_rows = []

        for sample in self.sample_names:
            result = self._analyze_sample(sample)
            if result is None:
                continue
            frag_parts.append(result['fragment_counts'])
            rt_parts.append(result['rt_dropoff'])
            len_parts.append(result['fragment_lengths'])
            if result.get('coverage') is not None:
                cov_parts.append(result['coverage'])
            summary_rows.append(result['summary'])

        self._fragment_counts = (pd.concat(frag_parts, ignore_index=True)
                                 if frag_parts else pd.DataFrame())
        self._rt_dropoff = (pd.concat(rt_parts, ignore_index=True)
                            if rt_parts else pd.DataFrame())
        self._fragment_lengths = (pd.concat(len_parts, ignore_index=True)
                                  if len_parts else pd.DataFrame())
        self._coverage = (pd.concat(cov_parts, ignore_index=True)
                          if cov_parts else pd.DataFrame())
        self._summary = (pd.DataFrame(summary_rows)
                         if summary_rows else pd.DataFrame())

    def _read_sample_csv(self, sample_name):
        """Read a per-sample stats CSV with chunked streaming for memory safety."""
        path = self.stats_dir / f'{sample_name}_stats.csv.bz2'
        if not path.exists():
            return None

        chunks = []
        # Only read columns we need; stream in 500k-row chunks
        available_cols = pd.read_csv(path, nrows=0).columns.tolist()
        usecols = [c for c in _USECOLS if c in available_cols]

        for chunk in pd.read_csv(path, usecols=usecols, chunksize=500_000):
            chunks.append(chunk)

        if not chunks:
            return None
        return pd.concat(chunks, ignore_index=True)

    def _analyze_sample(self, sample_name):
        """Analyze a single sample. Returns dict of DataFrames or None."""
        df = self._read_sample_csv(sample_name)
        if df is None or df.empty:
            return None

        # Ensure boolean columns
        for col in ('5p_cover', '3p_cover'):
            if col in df.columns:
                df[col] = df[col].astype(bool)

        # Classify fragments
        df['fragment_type'] = df.apply(
            lambda r: _FRAG_TYPES.get(
                (bool(r['5p_cover']), bool(r['3p_cover'])), 'degraded'),
            axis=1
        )

        count = self.count_col

        # ---- 1. Fragment counts per tRNA ----
        # Pivot: rows=tRNA, columns=fragment_type, values=sum(count)
        frag_pivot = (df.groupby(['tRNA_annotation', 'fragment_type'])[count]
                      .sum().unstack(fill_value=0).reset_index())

        # Ensure all 4 type columns exist
        for ft in ('full_length', 'rt_dropoff', '5p_tRF', 'degraded'):
            if ft not in frag_pivot.columns:
                frag_pivot[ft] = 0

        frag_pivot['total_reads'] = (frag_pivot['full_length'] +
                                     frag_pivot['rt_dropoff'] +
                                     frag_pivot['5p_tRF'] +
                                     frag_pivot['degraded'])

        # Filter low-count tRNAs
        frag_pivot = frag_pivot[frag_pivot['total_reads'] >= self.min_reads].copy()
        if frag_pivot.empty:
            return None

        # Add tRNA metadata from first occurrence
        meta_cols = ['tRNA_annotation', 'amino_acid', 'codon', 'anticodon']
        meta_cols = [c for c in meta_cols if c in df.columns]
        meta = df[meta_cols].drop_duplicates(subset=['tRNA_annotation'])
        frag_pivot = frag_pivot.merge(meta, on='tRNA_annotation', how='left')

        # Fractions
        total = frag_pivot['total_reads']
        for ft in ('full_length', 'rt_dropoff', '5p_tRF', 'degraded'):
            frag_pivot[f'frac_{ft}'] = frag_pivot[ft] / total

        # Integrity score
        denom = frag_pivot['full_length'] + frag_pivot['rt_dropoff']
        frag_pivot['integrity_score'] = np.where(
            denom > 0, frag_pivot['full_length'] / denom, np.nan)

        frag_pivot['sample_name_unique'] = sample_name

        # Reorder columns
        col_order = [
            'sample_name_unique', 'tRNA_annotation',
            'amino_acid', 'codon', 'anticodon',
            'total_reads', 'full_length', 'rt_dropoff', '5p_tRF', 'degraded',
            'frac_full_length', 'frac_rt_dropoff', 'frac_5p_tRF', 'frac_degraded',
            'integrity_score',
        ]
        col_order = [c for c in col_order if c in frag_pivot.columns]
        frag_df = frag_pivot[col_order]

        # ---- 2. RT drop-off positions ----
        rt_reads = df[df['fragment_type'] == 'rt_dropoff'].copy()
        if not rt_reads.empty and 'align_5p_idx' in rt_reads.columns:
            rt_pos = (rt_reads.groupby(
                ['tRNA_annotation', 'tRNA_annotation_len', 'align_5p_idx']
            )[count].sum().reset_index().rename(columns={
                'align_5p_idx': 'position', count: 'rt_stop_count'
            }))
            # Per-tRNA fraction
            trna_totals = rt_pos.groupby('tRNA_annotation')['rt_stop_count'].transform('sum')
            rt_pos['rt_stop_fraction'] = rt_pos['rt_stop_count'] / trna_totals
            rt_pos['sample_name_unique'] = sample_name
        else:
            rt_pos = pd.DataFrame()

        # ---- 3. Fragment lengths ----
        len_rows = []

        # RT drop-off: length = tRNA_annotation_len - align_5p_idx + 1
        if not rt_reads.empty and 'align_5p_idx' in rt_reads.columns:
            rt_reads = rt_reads.copy()
            rt_reads['fragment_length'] = (
                rt_reads['tRNA_annotation_len'] - rt_reads['align_5p_idx'] + 1
            )
            rl = (rt_reads.groupby(
                ['tRNA_annotation', 'fragment_type', 'fragment_length']
            )[count].sum().reset_index().rename(columns={count: 'read_count'}))
            len_rows.append(rl)

        # 5' tRF: length = align_3p_idx
        trf_reads = df[df['fragment_type'] == '5p_tRF'].copy()
        if not trf_reads.empty and 'align_3p_idx' in trf_reads.columns:
            trf_reads['fragment_length'] = trf_reads['align_3p_idx']
            rl = (trf_reads.groupby(
                ['tRNA_annotation', 'fragment_type', 'fragment_length']
            )[count].sum().reset_index().rename(columns={count: 'read_count'}))
            len_rows.append(rl)

        # Full-length: length = tRNA_annotation_len
        fl_reads = df[df['fragment_type'] == 'full_length'].copy()
        if not fl_reads.empty:
            fl_reads['fragment_length'] = fl_reads['tRNA_annotation_len']
            rl = (fl_reads.groupby(
                ['tRNA_annotation', 'fragment_type', 'fragment_length']
            )[count].sum().reset_index().rename(columns={count: 'read_count'}))
            len_rows.append(rl)

        if len_rows:
            frag_len_df = pd.concat(len_rows, ignore_index=True)
            # Per-tRNA+type fraction
            type_totals = frag_len_df.groupby(
                ['tRNA_annotation', 'fragment_type'])['read_count'].transform('sum')
            frag_len_df['fraction'] = frag_len_df['read_count'] / type_totals
            frag_len_df['sample_name_unique'] = sample_name
        else:
            frag_len_df = pd.DataFrame()

        # ---- 4. Coverage data (for Behrens/needle plots) ----
        cov_df = None
        if 'align_5p_idx' in df.columns and 'amino_acid' in df.columns:
            cov = (df.groupby(['amino_acid', 'align_5p_idx'])[count]
                   .sum().reset_index()
                   .rename(columns={'align_5p_idx': 'position', count: 'count'}))
            aa_len = (df.groupby('amino_acid')['tRNA_annotation_len']
                      .max().reset_index()
                      .rename(columns={'tRNA_annotation_len': 'max_tRNA_len'}))
            cov = cov.merge(aa_len, on='amino_acid', how='left')
            cov['sample_name_unique'] = sample_name
            cov_df = cov

        # ---- 5. Sample summary ----
        n_total = int(frag_df['total_reads'].sum())
        n_fl = int(frag_df['full_length'].sum())
        n_rt = int(frag_df['rt_dropoff'].sum())
        n_5p = int(frag_df['5p_tRF'].sum())
        n_deg = int(frag_df['degraded'].sum())

        summary_row = {
            'sample_name_unique': sample_name,
            'N_total_aligned': n_total,
            'N_full_length': n_fl,
            'N_rt_dropoff': n_rt,
            'N_5p_tRF': n_5p,
            'N_degraded': n_deg,
            'pct_full_length': round(n_fl / n_total * 100, 2) if n_total else 0,
            'pct_rt_dropoff': round(n_rt / n_total * 100, 2) if n_total else 0,
            'pct_5p_tRF': round(n_5p / n_total * 100, 2) if n_total else 0,
            'pct_degraded': round(n_deg / n_total * 100, 2) if n_total else 0,
            'median_integrity': round(float(frag_df['integrity_score'].median()), 4)
                if frag_df['integrity_score'].notna().any() else np.nan,
            'n_flagged_high_dropoff': int(
                (frag_df['frac_rt_dropoff'] > 0.5).sum()),
        }

        return {
            'fragment_counts': frag_df,
            'rt_dropoff': rt_pos,
            'fragment_lengths': frag_len_df,
            'coverage': cov_df,
            'summary': summary_row,
        }

    def get_fragment_counts(self, level='transcript'):
        """Return fragment counts aggregated at the given level.

        Parameters
        ----------
        level : str
            'transcript' (per-tRNA), 'codon', or 'aa'.
        """
        if self._fragment_counts is None or self._fragment_counts.empty:
            return pd.DataFrame()

        if level == 'transcript':
            return self._fragment_counts

        group_col = {'codon': 'codon', 'aa': 'amino_acid'}.get(level)
        if group_col is None or group_col not in self._fragment_counts.columns:
            return self._fragment_counts

        num_cols = ['total_reads', 'full_length', 'rt_dropoff', '5p_tRF', 'degraded']
        agg = (self._fragment_counts
               .groupby(['sample_name_unique', group_col])[num_cols]
               .sum().reset_index())

        total = agg['total_reads']
        for ft in ('full_length', 'rt_dropoff', '5p_tRF', 'degraded'):
            agg[f'frac_{ft}'] = agg[ft] / total

        denom = agg['full_length'] + agg['rt_dropoff']
        agg['integrity_score'] = np.where(denom > 0, agg['full_length'] / denom, np.nan)

        return agg

    def flag_unusual_fragments(self, rt_dropoff_threshold=0.5):
        """Return tRNAs with unusually high RT drop-off fraction."""
        if self._fragment_counts is None or self._fragment_counts.empty:
            return pd.DataFrame()
        mask = self._fragment_counts['frac_rt_dropoff'] > rt_dropoff_threshold
        return self._fragment_counts[mask].copy()

    def export(self, output_dir, write_csv=False):
        """Write results to output directory.

        Default format is Parquet. CSV copies written if ``write_csv=True``.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        datasets = {
            'fragment_counts': self._fragment_counts,
            'rt_dropoff_positions': self._rt_dropoff,
            'fragment_lengths': self._fragment_lengths,
            'fragment_summary': self._summary,
            'coverage_data': self._coverage,
        }

        for name, df in datasets.items():
            if df is None or df.empty:
                continue
            _save_df(df, output_dir / name, write_csv)


def _save_df(df, path_stem, write_csv=False):
    """Save DataFrame as parquet (and optionally CSV)."""
    try:
        df.to_parquet(f'{path_stem}.parquet', index=False)
    except Exception:
        df.to_csv(f'{path_stem}.csv', index=False)
        return
    if write_csv:
        df.to_csv(f'{path_stem}.csv', index=False)
