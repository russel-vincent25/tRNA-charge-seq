"""
Modification QC Report Generator
=================================

Generates modification_report.html from modification analysis outputs.
Self-contained Plotly HTML dashboard — follows the same pattern as
``trnaseq.qc.report.QCReportGenerator``.

Panels (in order):
    1. Per-sample modification counts (stacked bar: consensus vs non-consensus)
    2. Mismatch signature distribution (grouped bar by condition, using
       ``dominant_pattern`` column e.g. "T->C", "G->T"; falls back to
       ``modification`` column when ``dominant_pattern`` is absent)
    3. Modification landscape heatmap (tRNA x absolute position mapped to longest tRNA)
    4. Modification call reproducibility (stacked bar grouped by condition)
    5. Synthetic control error rate (only if synthetics present)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

import plotly.graph_objects as go
import plotly.offline
from plotly.subplots import make_subplots

from trnaseq.qc._common import (
    ReportContext,
    fig_to_div,
    render_html_shell,
    render_panel,
)


class ModificationReportGenerator:
    """Generate interactive HTML dashboard for modification QC."""

    MAX_PAIRS = 30
    MAX_TRNAS_HEATMAP = 50
    N_POS_BINS = 20  # kept for reference; absolute positions are used in panels 3 & 5

    def __init__(
        self,
        per_sample_calls: Dict[str, pd.DataFrame],
        aggregated_calls: Optional[pd.DataFrame] = None,
        consensus_calls: Optional[pd.DataFrame] = None,
        replicate_groups: Optional[Dict[str, List[str]]] = None,
        ref_dict: Optional[Dict[str, dict]] = None,
        summary_df: Optional[pd.DataFrame] = None,
        source_prefixes: Optional[dict] = None,
        context: Optional[ReportContext] = None,
    ):
        """
        Parameters
        ----------
        per_sample_calls : dict
            {sample_name_unique: calls DataFrame}.
        aggregated_calls : DataFrame, optional
            Output of ReplicateAggregator.aggregate() (all sites >=1 rep).
        consensus_calls : DataFrame, optional
            Subset where consensus_call==True.
        replicate_groups : dict, optional
            {condition: [sample_name_unique, ...]}.
        ref_dict : dict, optional
            {trna_name: {'seq': str, 'seq_len': int}} for 3' alignment.
        summary_df : DataFrame, optional
            modification_summary (one row per sample).
        source_prefixes : dict, optional
            {prefix: category} for classifying tRNA sources (e.g. synthetic).
        """
        self.per_sample_calls = per_sample_calls
        self.aggregated_calls = aggregated_calls
        self.consensus_calls = consensus_calls
        self.replicate_groups = replicate_groups
        self.ref_dict = ref_dict or {}
        self.summary_df = summary_df
        self.source_prefixes = source_prefixes
        self.context = context

        # Build position-mapping lookup once; reused by panels 3 & 5.
        self._pos_map_cache: Optional[Dict[str, np.ndarray]] = None
        self._max_len_cache: Optional[int] = None

    # ------------------------------------------------------------------
    # Main entry
    # ------------------------------------------------------------------

    def generate_html_report(self, output_path):
        """Generate self-contained HTML report at *output_path*."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        panels = []

        # 1. Per-sample counts (stacked bar)
        p1 = self._panel_per_sample_counts()
        if p1:
            panels.append(p1)

        # 2. Modification type distribution
        p2 = self._panel_modification_type_distribution()
        if p2:
            panels.append(p2)

        # 3. Modification landscape heatmap (mismatch rate)
        p3 = self._panel_modification_landscape()
        if p3:
            panels.append(p3)

        # 4. Modification call reproducibility (grouped by condition)
        p4 = self._panel_call_reproducibility()
        if p4:
            panels.append(p4)

        # 5. Synthetic control error rate (only if synthetics present)
        p5 = self._panel_synthetic_false_positives()
        if p5:
            panels.append(p5)

        body = '\n'.join(panels) if panels else '<p>No modification data available for report.</p>'

        html = render_html_shell(
            title='Modification QC Report',
            h1_text='Modification Analysis QC Dashboard',
            accent_color='#6c5ce7',
            body=body,
            context=self.context,
        )

        output_path.write_text(html)
        return str(output_path)

    # ---- helpers ----

    def _fig_to_div(self, fig):
        return fig_to_div(fig)

    @staticmethod
    def _discrete_palette() -> List[str]:
        """Return the standard discrete color palette used across all panels."""
        return [
            '#6c5ce7', '#0984e3', '#00b894', '#fdcb6e', '#e17055',
            '#fd79a8', '#a29bfe', '#55efc4', '#ffeaa7', '#fab1a0',
            '#74b9ff', '#81ecec', '#dfe6e9', '#636e72', '#2d3436',
        ]

    @staticmethod
    def _isoacceptor_sort_key(trna_name: str) -> str:
        """Sort key: amino acid + anticodon extracted from FASTA name."""
        parts = trna_name.split('-')
        aa = parts[1] if len(parts) > 1 else ''
        ac = parts[2] if len(parts) > 2 else ''
        return f'{aa}-{ac}'

    def _normalize_position(self, row, pos_col='position') -> Optional[float]:
        """Return position as percent of tRNA length (0-100), or None if unknown."""
        info = self.ref_dict.get(row['trna_name'])
        seq_len = info['seq_len'] if info else None
        if seq_len and seq_len > 0:
            return float(row[pos_col]) / seq_len * 100.0
        return None

    @staticmethod
    def _bin_normalized_pos(norm_pos: float, n_bins: int = 20) -> int:
        """Map 0-100% position to bin index [0, n_bins-1]."""
        return min(int(norm_pos / 100.0 * n_bins), n_bins - 1)

    @staticmethod
    def _bin_labels(n_bins: int = 20) -> List[str]:
        step = 100 // n_bins
        return [f'{i * step}-{(i + 1) * step}%' for i in range(n_bins)]

    # ------------------------------------------------------------------
    # Position mapping: absolute coords onto longest tRNA
    # ------------------------------------------------------------------

    def _build_pos_map(self, df: pd.DataFrame):
        """Build (and cache) per-tRNA position -> longest-tRNA-index lookup.

        Uses the nearest-percentile algorithm from TRNA_plot coverage:
            for each tRNA of length L, map its positions [0..L-1] to the
            nearest integer in [0..max_len-1] using np.percentile with
            np.linspace(0, 100, L) as quantiles over np.arange(max_len).

        If ref_dict is available, tRNA lengths come from 'seq_len'.
        Otherwise lengths are estimated as ceil(max_position_in_data * 1.5)
        per tRNA (modifications are typically in the first 2/3 of the tRNA).

        Returns
        -------
        pos_map : dict  {trna_name: ndarray of shape (trna_len,)}
            pos_map[t][p] gives the mapped index in [0, max_len-1].
        max_len : int
        """
        if self._pos_map_cache is not None:
            return self._pos_map_cache, self._max_len_cache

        # Gather tRNA lengths
        trna_names = df['trna_name'].unique() if 'trna_name' in df.columns else []
        lengths: Dict[str, int] = {}
        for t in trna_names:
            if self.ref_dict and t in self.ref_dict:
                info = self.ref_dict[t]
                l = info.get('seq_len') or (len(info['seq']) if info.get('seq') else None)
                if l and l > 0:
                    lengths[t] = int(l)
            if t not in lengths:
                # Estimate: max position seen * 1.5 (positions are usually in first 2/3)
                max_p = df[df['trna_name'] == t]['position'].max() if not df.empty else 0
                if pd.notna(max_p) and max_p > 0:
                    lengths[t] = max(1, int(np.ceil(max_p * 1.5)))
                else:
                    lengths[t] = 1

        if not lengths:
            self._pos_map_cache = {}
            self._max_len_cache = 1
            return self._pos_map_cache, self._max_len_cache

        max_len = max(lengths.values())
        max_len = max(max_len, 1)

        pos_map: Dict[str, np.ndarray] = {}
        ref_indices = np.arange(max_len)

        for t, tlen in lengths.items():
            if tlen == max_len:
                # Identity mapping
                pos_map[t] = ref_indices.copy()
            elif tlen == 1:
                pos_map[t] = np.array([0], dtype=int)
            else:
                quantiles = np.linspace(0, 100, tlen)
                mapped = np.percentile(ref_indices, quantiles, method='nearest').astype(int)
                pos_map[t] = mapped

        self._pos_map_cache = pos_map
        self._max_len_cache = max_len
        return pos_map, max_len

    def _map_position_to_longest(self, trna_name: str, position, pos_map: dict, max_len: int) -> Optional[int]:
        """Map a tRNA position to the longest-tRNA coordinate.

        Parameters
        ----------
        trna_name : str
        position : int or float — 0-based position within the tRNA
        pos_map : dict — output of _build_pos_map
        max_len : int — output of _build_pos_map

        Returns
        -------
        Mapped integer index in [0, max_len-1], or None if unmappable.
        """
        mapping = pos_map.get(trna_name)
        if mapping is None:
            return None
        pos = int(position)
        if 0 <= pos < len(mapping):
            return int(mapping[pos])
        return None

    def _apply_position_mapping(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add '_mapped_pos' column to df using absolute position mapping.

        Rows where mapping fails are dropped.
        """
        df = df.copy()
        pos_map, max_len = self._build_pos_map(df)

        mapped = []
        for _, row in df.iterrows():
            mp = self._map_position_to_longest(
                row['trna_name'], row['position'], pos_map, max_len
            )
            mapped.append(mp)

        df['_mapped_pos'] = mapped
        df = df.dropna(subset=['_mapped_pos'])
        df['_mapped_pos'] = df['_mapped_pos'].astype(int)
        return df

    # ------------------------------------------------------------------
    # Panel 1: per-sample modification counts (stacked bar)
    # ------------------------------------------------------------------

    def _panel_per_sample_counts(self):
        if not self.per_sample_calls:
            return ''

        samples = sorted(self.per_sample_calls.keys())
        total_counts = []

        for s in samples:
            df = self.per_sample_calls[s]
            total_counts.append(len(df) if df is not None and not df.empty else 0)

        # Determine consensus counts per sample
        has_consensus = (
            self.summary_df is not None
            and 'consensus_calls' in self.summary_df.columns
        )

        consensus_counts = []
        if has_consensus:
            for s in samples:
                row = self.summary_df[self.summary_df['sample_name_unique'] == s]
                consensus_counts.append(
                    int(row['consensus_calls'].iloc[0]) if not row.empty else 0
                )
        else:
            consensus_counts = [0] * len(samples)

        non_consensus_counts = [
            max(0, t - c) for t, c in zip(total_counts, consensus_counts)
        ]

        fig = go.Figure()

        if has_consensus:
            # Stacked: bottom = consensus (green), top = non-consensus (orange)
            fig.add_trace(go.Bar(
                x=samples,
                y=consensus_counts,
                name='Consensus calls',
                marker_color='#00b894',
            ))
            fig.add_trace(go.Bar(
                x=samples,
                y=non_consensus_counts,
                name='Non-consensus calls',
                marker_color='#fdcb6e',
            ))
            barmode = 'stack'
        else:
            # No consensus info — show total only
            fig.add_trace(go.Bar(
                x=samples,
                y=total_counts,
                name='Total calls',
                marker_color='#0984e3',
            ))
            barmode = 'group'

        fig.update_layout(
            title='Modification Calls per Sample',
            yaxis_title='Number of modifications',
            barmode=barmode,
            height=max(350, len(samples) * 20 + 200),
        )
        div = self._fig_to_div(fig)
        return render_panel('Per-Sample Modification Counts', div, anchor='mod-counts')

    # ------------------------------------------------------------------
    # Panel 2: mismatch signature distribution — grouped by condition
    # ------------------------------------------------------------------

    def _panel_modification_type_distribution(self):
        """Grouped bar chart of dominant mismatch pattern counts, one bar group per condition.

        Uses the ``dominant_pattern`` column (e.g. "T->C", "G->T", "A->T")
        when available, because most modification calls are not matched to a
        MODOMICS entry and the raw mismatch signature is the most informative
        signal.  Falls back to the ``modification`` column when
        ``dominant_pattern`` is absent (legacy behaviour).

        Conditions come from self.replicate_groups keys.  Each condition's
        counts are summed from the per-sample calls of its member samples.
        Falls back to a single bar group ('All samples') when
        replicate_groups is not available.

        Patterns are sorted by total count (descending) so the most common
        signatures appear first.
        """
        if not self.per_sample_calls:
            return ''

        # Decide which column to use: prefer dominant_pattern, fall back to modification
        use_col = None
        for sample, df in self.per_sample_calls.items():
            if df is None or df.empty:
                continue
            if 'dominant_pattern' in df.columns:
                use_col = 'dominant_pattern'
            elif 'modification' in df.columns:
                use_col = 'modification'
            break  # column availability is the same for all samples

        if use_col is None:
            return ''

        using_patterns = (use_col == 'dominant_pattern')

        # Collect per-sample counts, dropping empty/NaN values
        all_labels: set = set()
        sample_label_counts: Dict[str, pd.Series] = {}

        for sample, df in self.per_sample_calls.items():
            if df is None or df.empty or use_col not in df.columns:
                continue
            series = df[use_col].dropna()
            series = series[series.astype(str).str.strip() != '']
            if series.empty:
                continue
            counts = series.value_counts()
            sample_label_counts[sample] = counts
            all_labels.update(counts.index.tolist())

        if not all_labels:
            return ''

        # Build condition -> {label: count} mapping
        if self.replicate_groups:
            conditions = list(self.replicate_groups.keys())
            cond_counts: Dict[str, Dict[str, int]] = {}
            for cond, members in self.replicate_groups.items():
                agg: Dict[str, int] = {lbl: 0 for lbl in all_labels}
                for sample in members:
                    counts = sample_label_counts.get(sample)
                    if counts is not None:
                        for lbl in all_labels:
                            agg[lbl] += int(counts.get(lbl, 0))
                cond_counts[cond] = agg
        else:
            conditions = ['All samples']
            total: Dict[str, int] = {lbl: 0 for lbl in all_labels}
            for counts in sample_label_counts.values():
                for lbl in all_labels:
                    total[lbl] += int(counts.get(lbl, 0))
            cond_counts = {'All samples': total}

        # Sort labels by total count descending
        label_totals = {
            lbl: sum(cond_counts[c].get(lbl, 0) for c in conditions)
            for lbl in all_labels
        }
        labels = sorted(all_labels, key=lambda x: -label_totals[x])

        # One color per condition from the discrete palette
        palette = self._discrete_palette()

        fig = go.Figure()
        for i, cond in enumerate(conditions):
            y_vals = [cond_counts[cond].get(lbl, 0) for lbl in labels]
            fig.add_trace(go.Bar(
                x=labels,
                y=y_vals,
                name=cond,
                marker_color=palette[i % len(palette)],
                text=y_vals,
                textposition='auto',
            ))

        if using_patterns:
            title = 'Mismatch Signature Distribution by Condition'
            xaxis_title = 'Mismatch pattern (reference \u2192 observed)'
            subtitle = (
                'Mismatch patterns indicate the reference\u2192observed substitution at '
                'modification sites. Common patterns: T\u2192C (pseudouridine/dihydrouridine), '
                'G\u2192T (m\u2077G), A\u2192T (i\u2076A/t\u2076A). '
                'Patterns are sorted by total count across all conditions.'
            )
        else:
            title = 'Modification Type Distribution by Condition'
            xaxis_title = 'Modification type'
            subtitle = (
                'Counts of each modification type per condition. '
                'Reveals whether certain types are specific to particular '
                'RT enzymes or conditions.'
            )

        fig.update_layout(
            title=title,
            xaxis_title=xaxis_title,
            yaxis_title='Number of sites',
            barmode='group',
            height=max(400, len(labels) * 18 + 200),
            legend_title='Condition',
        )
        div = self._fig_to_div(fig)

        heading = 'Mismatch Signature Distribution' if using_patterns else 'Modification Type Distribution'
        return render_panel(heading, div, anchor='mod-types', description=subtitle)

    # ------------------------------------------------------------------
    # Panel 3: modification landscape heatmap (mismatch rate, absolute positions)
    # ------------------------------------------------------------------

    def _panel_modification_landscape(self):
        """Heatmap of mean mismatch rate across tRNA positions x tRNA names.

        Rows = tRNA species (sorted by isoacceptor), columns = absolute
        position index mapped to the longest tRNA in the dataset (same
        nearest-percentile algorithm as the Behrens/Needle coverage plots).
        Cell color encodes mean mismatch rate.
        """
        df = self.aggregated_calls if (
            self.aggregated_calls is not None and not self.aggregated_calls.empty
        ) else self.consensus_calls

        if df is None or df.empty:
            return ''

        rate_col = None
        for candidate in ('mean_mismatch_rate', 'mismatch_rate'):
            if candidate in df.columns:
                rate_col = candidate
                break
        if rate_col is None or 'trna_name' not in df.columns or 'position' not in df.columns:
            return ''

        df = df.copy()
        df = self._apply_position_mapping(df)

        if df.empty:
            return ''

        _, max_len = self._build_pos_map(df)

        # Sort tRNAs by isoacceptor
        trnas = sorted(df['trna_name'].unique(), key=self._isoacceptor_sort_key)
        if len(trnas) > self.MAX_TRNAS_HEATMAP:
            trnas = trnas[:self.MAX_TRNAS_HEATMAP]
            df = df[df['trna_name'].isin(trnas)]

        # Aggregate: mean mismatch rate per (trna, mapped position)
        agg = df.groupby(['trna_name', '_mapped_pos'])[rate_col].mean().reset_index()

        matrix = np.full((len(trnas), max_len), np.nan)
        trna_idx = {t: i for i, t in enumerate(trnas)}

        for _, row in agg.iterrows():
            ti = trna_idx.get(row['trna_name'])
            pi = int(row['_mapped_pos'])
            if ti is not None and 0 <= pi < max_len:
                matrix[ti, pi] = row[rate_col]

        x_labels = list(range(max_len))

        fig = go.Figure(go.Heatmap(
            z=matrix,
            x=x_labels,
            y=trnas,
            colorscale='YlOrRd',
            zmin=0,
            zmax=0.5,
            colorbar_title='Mean mismatch rate',
            hoverongaps=False,
        ))
        fig.update_layout(
            title="Modification Landscape (Mean Mismatch Rate by Position)",
            xaxis_title="5' to 3' index (mapped to longest tRNA)",
            yaxis_title='tRNA (sorted by isoacceptor)',
            height=max(400, len(trnas) * 18 + 150),
        )
        div = self._fig_to_div(fig)
        landscape_desc = (
            'Mean mismatch rate per tRNA and absolute position. Positions from all '
            'tRNA lengths are mapped to the longest tRNA using nearest-percentile '
            'interpolation (consistent with Behrens/Needle coverage plots). '
            'Bright cells indicate positions with elevated substitution signal '
            '(likely RNA modifications). Each row is one tRNA species; '
            "x-axis spans 5′ to 3′ end."
        )
        return render_panel(
            'Modification Landscape', div, anchor='mod-landscape',
            description=landscape_desc,
        )

    # ------------------------------------------------------------------
    # Panel 4: modification call reproducibility (grouped by condition)
    # ------------------------------------------------------------------

    def _panel_call_reproducibility(self):
        """Stacked bar showing number of modification sites detected in
        1, 2, 3, ... replicates, with one bar group per condition.

        Uses ``self.replicate_groups`` ({condition: [sample_names]}) to assign
        each site to a condition based on sample membership.  Falls back to a
        single 'All samples' bar when replicate_groups is unavailable.

        Condition labels are the keys of self.replicate_groups, which are
        sample names with the -R\\d+ suffix stripped (e.g. "TGIRT-37", "SSIV-55").
        """
        if self.aggregated_calls is None or self.aggregated_calls.empty:
            return ''
        if 'n_replicates_detected' not in self.aggregated_calls.columns:
            return ''

        agg = self.aggregated_calls.copy()

        # Determine max replicates
        max_reps = int(agg['n_replicates_detected'].max())
        if max_reps < 1:
            return ''

        rep_levels = list(range(1, max_reps + 1))

        # Color gradient: 1 rep = yellow, 2 = orange, 3 = red, 4+ = darker
        level_colors = {
            1: '#fdcb6e',
            2: '#e17055',
            3: '#d63031',
            4: '#6c5ce7',
            5: '#2d3436',
        }

        if self.replicate_groups:
            # sample_name in aggregated_calls IS already the condition group
            if 'sample_name' in agg.columns:
                agg['_condition'] = agg['sample_name']
            else:
                agg['_condition'] = 'All samples'

            all_conds_in_data = set(agg['_condition'].unique())
            conditions = [c for c in self.replicate_groups.keys() if c in all_conds_in_data]
            for c in sorted(all_conds_in_data - set(conditions)):
                conditions.append(c)

            fig = go.Figure()
            for rep in rep_levels:
                y_vals = []
                for cond in conditions:
                    subset = agg[
                        (agg['_condition'] == cond) &
                        (agg['n_replicates_detected'] == rep)
                    ]
                    # Count unique (trna_name, position) sites
                    if 'trna_name' in agg.columns and 'position' in agg.columns:
                        n_sites = subset[['trna_name', 'position']].drop_duplicates().shape[0]
                    else:
                        n_sites = len(subset)
                    y_vals.append(n_sites)

                fig.add_trace(go.Bar(
                    x=conditions,
                    y=y_vals,
                    name=f'{rep} replicate{"s" if rep > 1 else ""}',
                    marker_color=level_colors.get(rep, '#b2bec3'),
                ))

            x_label = 'Condition'
        else:
            # Fallback: single bar group across all data
            conditions = ['All samples']
            fig = go.Figure()
            for rep in rep_levels:
                subset = agg[agg['n_replicates_detected'] == rep]
                if 'trna_name' in agg.columns and 'position' in agg.columns:
                    n_sites = subset[['trna_name', 'position']].drop_duplicates().shape[0]
                else:
                    n_sites = len(subset)
                fig.add_trace(go.Bar(
                    x=conditions,
                    y=[n_sites],
                    name=f'{rep} replicate{"s" if rep > 1 else ""}',
                    marker_color=level_colors.get(rep, '#b2bec3'),
                ))
            x_label = ''

        fig.update_layout(
            title='Modification Call Reproducibility by Condition',
            xaxis_title=x_label,
            yaxis_title='Number of sites',
            barmode='stack',
            height=max(400, 350),
        )
        div = self._fig_to_div(fig)
        return render_panel('Modification Call Reproducibility', div, anchor='mod-reproducibility')

    # ------------------------------------------------------------------
    # Panel 5: synthetic control error rate
    # ------------------------------------------------------------------

    def _panel_synthetic_false_positives(self):
        """Bar chart of mean mismatch rate on synthetic control tRNAs per sample.

        Returns '' (skipped) if no synthetic tRNAs are detected in any sample.
        """
        if not self.source_prefixes or not self.per_sample_calls:
            return ''

        from trnaseq.charge.quantifier import classify_trna_source

        samples = sorted(self.per_sample_calls.keys())
        syn_error_rates = {}  # {sample: mean mismatch rate on synthetics}
        any_synthetics = False

        rate_col = None  # determined from first non-empty df

        for sample in samples:
            calls_df = self.per_sample_calls[sample]
            if calls_df is None or calls_df.empty:
                syn_error_rates[sample] = np.nan
                continue

            trna_col = 'trna_name' if 'trna_name' in calls_df.columns else None
            if trna_col is None:
                syn_error_rates[sample] = np.nan
                continue

            # Detect rate column once
            if rate_col is None:
                for candidate in ('mismatch_rate', 'mean_mismatch_rate'):
                    if candidate in calls_df.columns:
                        rate_col = candidate
                        break

            calls_df = calls_df.copy()
            calls_df['_source'] = calls_df[trna_col].apply(
                lambda x: classify_trna_source(x, self.source_prefixes)
            )
            syn_calls = calls_df[calls_df['_source'] == 'synthetic']

            if syn_calls.empty or rate_col is None or rate_col not in syn_calls.columns:
                syn_error_rates[sample] = np.nan
            else:
                any_synthetics = True
                syn_error_rates[sample] = syn_calls[rate_col].mean()

        # Skip panel entirely if no synthetic tRNAs detected in any sample
        if not any_synthetics:
            return ''

        y_vals = [syn_error_rates.get(s, np.nan) for s in samples]

        fig = go.Figure(go.Bar(
            x=samples,
            y=y_vals,
            marker_color='#e17055',
            text=[f'{v:.4f}' if not np.isnan(v) else 'N/A' for v in y_vals],
            textposition='auto',
        ))

        # Background error rate reference line
        if self.summary_df is not None and 'background_error_rate' in self.summary_df.columns:
            bg_rate = float(self.summary_df['background_error_rate'].mean())
            fig.add_hline(
                y=bg_rate,
                line_dash='dash',
                line_color='#636e72',
                annotation_text=f'Background error rate: {bg_rate:.4f}',
                annotation_position='top right',
            )

        fig.update_layout(
            title='Error Rate on Synthetic Controls',
            xaxis_title='Sample',
            yaxis_title='Mean mismatch rate (synthetic tRNAs)',
            height=max(350, len(samples) * 20 + 200),
        )

        div = self._fig_to_div(fig)
        return render_panel('Synthetic Control Error Rate', div, anchor='mod-synthetic')
