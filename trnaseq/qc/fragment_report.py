"""
Fragment QC Report Generator
==============================

Generates fragment_report.html — an interactive Plotly dashboard
summarizing tRNA fragment analysis across samples.

Panels:
    1. Fragment type stacked bar
    2. Fragment composition by tRNA — stacked bar of fragment type fractions
       per tRNA for the top 30 most-covered tRNAs, with a dropdown to switch
       between conditions (each condition averaged across its replicates)
    3. Per-AA integrity box plots with condition dropdown
    4. Behrens coverage plot (if coverage data provided)
    5. Needle coverage plot (if coverage data provided)
    6. Synthetic tRNA integrity (if source_prefixes provided)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional

import plotly.graph_objects as go
import plotly.offline

from trnaseq.qc.report import _discrete_palette
from trnaseq.qc._common import (
    ReportContext,
    fig_to_div,
    render_html_shell,
    render_panel,
)


def _isoacceptor_sort_key(trna_name):
    """Sort key for tRNA names by amino acid then anticodon."""
    parts = trna_name.split('-')
    if len(parts) >= 3:
        return (parts[1], parts[2])
    return (trna_name, '')


# Standard amino acid order — reverse alphabetical (matches original TRNA_plot)
_AA_ORDER = [
    'Val', 'Tyr', 'Trp', 'Thr', 'Ser', 'Pro', 'Phe', 'Met',
    'Lys', 'Leu', 'Ile', 'His', 'Gly', 'Glu', 'Gln', 'Cys',
    'Asp', 'Asn', 'Arg', 'Ala',
]


def _aa_palette():
    """Return a {amino_acid: hex_color} mapping using tab20."""
    palette = _discrete_palette(20)
    return {aa: palette[i] for i, aa in enumerate(_AA_ORDER)}


class FragmentReportGenerator:
    """Generate interactive HTML dashboard for fragment analysis QC."""

    MAX_TRNAS_HEATMAP = 50

    def __init__(
        self,
        fragment_counts_df: pd.DataFrame,
        rt_dropoff_df: pd.DataFrame,
        fragment_lengths_df: pd.DataFrame,
        fragment_summary_df: pd.DataFrame,
        sample_df: pd.DataFrame,
        coverage_df: Optional[pd.DataFrame] = None,
        source_prefixes: Optional[dict] = None,
        context: Optional[ReportContext] = None,
    ):
        """
        Parameters
        ----------
        fragment_counts_df : DataFrame
            Per-tRNA per-sample fragment counts with integrity_score.
        rt_dropoff_df : DataFrame
            RT drop-off positional data (position, rt_stop_fraction).
        fragment_lengths_df : DataFrame
            Fragment length distributions.
        fragment_summary_df : DataFrame
            Per-sample summary (N_full_length, N_rt_dropoff, etc.).
        sample_df : DataFrame
            Sample metadata.
        coverage_df : DataFrame, optional
            Per-sample per-AA per-position coverage data from FragmentAnalyser.
            Columns: sample_name_unique, amino_acid, position, count, max_tRNA_len.
        source_prefixes : dict, optional
            {prefix: category} for classifying tRNA sources (e.g. synthetic).
        """
        self.frag_counts = fragment_counts_df.copy()
        self.rt_dropoff = rt_dropoff_df.copy()
        self.frag_lengths = fragment_lengths_df.copy()
        self.frag_summary = fragment_summary_df.copy()
        self.sample_df = sample_df.copy()
        self.coverage_df = coverage_df.copy() if coverage_df is not None else None
        self.source_prefixes = source_prefixes
        self.context = context

    def generate_html_report(self, output_path) -> str:
        """Generate self-contained HTML report."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        panels = []

        p1 = self._panel_fragment_type_bar()
        if p1:
            panels.append(p1)

        p2 = self._panel_fragment_composition_by_trna()
        if p2:
            panels.append(p2)

        p3 = self._panel_integrity_violin()
        if p3:
            panels.append(p3)

        p4 = self._panel_behrens_coverage()
        if p4:
            panels.append(p4)

        p5 = self._panel_needle_coverage()
        if p5:
            panels.append(p5)

        p6 = self._panel_synthetic_integrity()
        if p6:
            panels.append(p6)

        body = '\n'.join(panels)

        html = render_html_shell(
            title='tRNA Fragment Report',
            h1_text='tRNA Fragment Dashboard',
            accent_color='#0984e3',
            body=body,
            context=self.context,
        )

        output_path.write_text(html)
        return str(output_path)

    def _fig_to_div(self, fig):
        return fig_to_div(fig)

    # ---- panel 1: fragment type stacked bar ----

    def _panel_fragment_type_bar(self):
        """Stacked bar of fragment types per sample."""
        df = self.frag_summary
        if df.empty:
            return ''

        type_cols = {
            'N_full_length': ('Full-length', '#00b894'),
            'N_rt_dropoff': ('RT drop-off', '#0984e3'),
            'N_5p_tRF': ("5' tRF", '#fdcb6e'),
            'N_degraded': ('Degraded', '#d63031'),
        }

        present = {k: v for k, v in type_cols.items() if k in df.columns}
        if not present:
            return ''

        fig = go.Figure()
        for col, (label, color) in present.items():
            fig.add_trace(go.Bar(
                x=df['sample_name_unique'],
                y=df[col],
                name=label,
                marker_color=color,
            ))
        fig.update_layout(
            barmode='stack',
            title='Fragment Type Distribution per Sample',
            yaxis_title='Read count',
            height=max(400, len(df) * 25),
        )
        div = self._fig_to_div(fig)
        return render_panel('Fragment Types', div, anchor='frag-types')

    # ---- panel 2: fragment composition by tRNA ----

    TOP_TRNAS_COMPOSITION = 30

    def _panel_fragment_composition_by_trna(self):
        """Stacked bar of fragment type fractions per tRNA, with condition dropdown.

        Uses frac_* columns from frag_counts. Selects the top 30 tRNAs by
        total_reads across all samples, sorted by isoacceptor. For each
        condition (derived by stripping the trailing -R\\d+ suffix), the mean
        fraction per fragment type per tRNA is computed across replicates. A
        dropdown lets the user switch between conditions; the first condition is
        visible by default.
        """
        import re

        df = self.frag_counts
        if df.empty:
            return ''

        trna_col = 'tRNA_annotation' if 'tRNA_annotation' in df.columns else 'tRNA_name'
        frac_cols = ['frac_full_length', 'frac_rt_dropoff', 'frac_5p_tRF', 'frac_degraded']
        required = {trna_col, 'sample_name_unique', 'total_reads'} | set(frac_cols)
        if not required.issubset(df.columns):
            return ''

        df = df.copy()
        df['_condition'] = df['sample_name_unique'].apply(
            lambda s: re.sub(r'-R\d+$', '', str(s))
        )

        # Select top 30 tRNAs by total_reads summed across all samples
        total_by_trna = df.groupby(trna_col)['total_reads'].sum()
        top_trnas = total_by_trna.nlargest(self.TOP_TRNAS_COMPOSITION).index.tolist()
        # Sort by isoacceptor (AA then anticodon)
        top_trnas = sorted(top_trnas, key=_isoacceptor_sort_key)

        sub = df[df[trna_col].isin(top_trnas)]
        if sub.empty:
            return ''

        # Short display names: strip prefix, keep AA-anticodon-copy (everything after first '-')
        def _short_name(name):
            idx = name.find('-')
            return name[idx + 1:] if idx != -1 else name

        x_labels = [_short_name(t) for t in top_trnas]

        # Fragment type config: internal key → (display label, color)
        type_colors = {
            'frac_full_length': ('Full-length', '#00b894'),
            'frac_rt_dropoff':  ('RT drop-off',  '#0984e3'),
            'frac_5p_tRF':      ("5' tRF",        '#fdcb6e'),
            'frac_degraded':    ('Degraded',       '#d63031'),
        }

        conditions = sorted(sub['_condition'].unique())
        if not conditions:
            return ''

        fig = go.Figure()
        traces_per_condition = []

        for cond_idx, cond in enumerate(conditions):
            visible = (cond_idx == 0)
            cond_df = sub[sub['_condition'] == cond]

            # Mean fraction per fragment type per tRNA across replicates
            mean_fracs = (
                cond_df.groupby(trna_col)[frac_cols]
                .mean()
            )

            for ftype, (label, color) in type_colors.items():
                y_vals = [
                    float(mean_fracs.loc[t, ftype]) if t in mean_fracs.index else float('nan')
                    for t in top_trnas
                ]
                fig.add_trace(go.Bar(
                    x=x_labels,
                    y=y_vals,
                    name=label,
                    marker_color=color,
                    visible=visible,
                    showlegend=(cond_idx == 0),
                    hovertemplate=(
                        f'{label}<br>tRNA: %{{x}}<br>'
                        'Fraction: %{y:.3f}<extra></extra>'
                    ),
                ))
            traces_per_condition.append(4)

        # Build dropdown buttons
        total = sum(traces_per_condition)
        buttons = []
        offset = 0
        for cond_idx, cond in enumerate(conditions):
            n_t = traces_per_condition[cond_idx]
            vis = [False] * total
            for j in range(n_t):
                vis[offset + j] = True
            buttons.append(dict(
                label=cond,
                method='update',
                args=[
                    {'visible': vis},
                    {'title': f'Fragment Composition by tRNA \u2014 {cond}'},
                ],
            ))
            offset += n_t

        first_cond = conditions[0]
        fig.update_layout(
            barmode='stack',
            title=f'Fragment Composition by tRNA \u2014 {first_cond}',
            xaxis_title='tRNA',
            yaxis_title='Fraction',
            yaxis_range=[0, 1],
            height=550,
            updatemenus=[dict(
                type='dropdown',
                direction='down',
                x=1.0, xanchor='right',
                y=1.15, yanchor='top',
                buttons=buttons,
                active=0,
            )],
        )
        div = self._fig_to_div(fig)
        return render_panel('Fragment Composition by tRNA', div, anchor='frag-composition')

    # ---- panel 3: per-AA integrity violin ----

    def _panel_integrity_violin(self):
        """Box plots of integrity_score by amino acid, with a dropdown to select condition."""
        import re

        df = self.frag_counts
        if df.empty or 'integrity_score' not in df.columns:
            return ''

        if 'amino_acid' not in df.columns or 'sample_name_unique' not in df.columns:
            return ''

        aa_colors = _aa_palette()

        # Extract condition by stripping trailing -R\d+ replicate suffix
        def _condition(name):
            return re.sub(r'-R\d+$', '', str(name))

        df = df.copy()
        df['_condition'] = df['sample_name_unique'].apply(_condition)
        conditions = sorted(df['_condition'].unique())
        if not conditions:
            return ''

        # Collect all amino acids across the whole dataset (for consistent x-axis)
        all_aas = sorted(df['amino_acid'].dropna().unique())
        if not all_aas:
            return ''

        fig = go.Figure()
        traces_per_condition = []

        for c_idx, cond in enumerate(conditions):
            visible = (c_idx == 0)
            cond_df = df[df['_condition'] == cond]
            n_traces = 0
            for aa in all_aas:
                sub = cond_df[cond_df['amino_acid'] == aa]['integrity_score'].dropna()
                if sub.empty:
                    continue
                fig.add_trace(go.Box(
                    y=sub,
                    name=aa,
                    marker_color=aa_colors.get(aa, '#b2bec3'),
                    boxpoints='outliers',
                    visible=visible,
                    showlegend=False,
                    hovertemplate=f'{aa} | {cond}<br>Integrity: %{{y:.3f}}<extra></extra>',
                ))
                n_traces += 1
            traces_per_condition.append(n_traces)

        if not any(n > 0 for n in traces_per_condition):
            return ''

        # Build dropdown buttons
        total = sum(traces_per_condition)
        final_buttons = []
        offset = 0
        for c_idx, cond in enumerate(conditions):
            n_t = traces_per_condition[c_idx]
            vis = [False] * total
            for j in range(n_t):
                vis[offset + j] = True
            final_buttons.append(dict(
                label=cond,
                method='update',
                args=[{'visible': vis}, {'title': f'tRNA Integrity by Amino Acid — {cond}'}],
            ))
            offset += n_t

        first_cond = conditions[0]
        fig.update_layout(
            title=f'tRNA Integrity by Amino Acid — {first_cond}',
            yaxis_title='Integrity score',
            yaxis_range=[-0.05, 1.05],
            height=450,
            showlegend=False,
            updatemenus=[dict(
                type='dropdown',
                direction='down',
                x=1.0, xanchor='right',
                y=1.15, yanchor='top',
                buttons=final_buttons,
                active=0,
            )],
        )
        div = self._fig_to_div(fig)
        return render_panel('Integrity by Amino Acid', div, anchor='frag-integrity')

    # ---- coverage normalization helper ----

    def _build_coverage_matrices(self, sample_cov, aa_norm=False):
        """Build coverage matrices matching the original TRNA_plot algorithm.

        1. Map variable-length tRNAs to the global max length using
           nearest-percentile indexing.
        2. Place read-start counts at mapped 5' positions.
        3. Cumulate 5'→3' (left-to-right) to get coverage.
        4. Column-wise cumulate (stack) for Behrens-style plotting.

        Parameters
        ----------
        sample_cov : DataFrame
            Subset of coverage_df for a single sample.
            Columns: amino_acid, position, count, max_tRNA_len.
        aa_norm : bool
            If True, weight each amino acid equally at the 3' end.

        Returns
        -------
        cov_count : ndarray (n_aa, max_len)
            Per-AA coverage (cumulated 5'→3').
        cov_count_sum : ndarray (n_aa+1, max_len)
            Column-wise cumulated coverage (row 0 = zeros baseline).
        aa_ordered : list of str
            Amino acid names in plot order (reverse alphabetical).
        max_len : int
            Length of the x-axis.
        """
        if sample_cov.empty:
            return None, None, [], 0

        # Global max tRNA length across all AAs
        global_max_len = int(sample_cov['max_tRNA_len'].max())
        if global_max_len < 2:
            return None, None, [], 0

        # Build nearest-percentile length maps for each observed length
        observed_lens = sample_cov['max_tRNA_len'].unique()
        len_map = {}
        for tlen in observed_lens:
            tlen = int(tlen)
            if tlen not in len_map:
                len_map[tlen] = np.percentile(
                    np.arange(global_max_len),
                    np.linspace(0, 100, tlen),
                    method='nearest',
                ).astype(int)

        # Order amino acids: use _AA_ORDER for known, then extras
        all_aas = set(sample_cov['amino_acid'].unique())
        aa_ordered = [aa for aa in _AA_ORDER if aa in all_aas]
        for extra in sorted(all_aas - set(_AA_ORDER)):
            aa_ordered.append(extra)
        aa_index = {aa: i for i, aa in enumerate(aa_ordered)}

        # Build count matrix: place read-start counts at mapped positions
        n_aa = len(aa_ordered)
        cov_count = np.zeros((n_aa, global_max_len))

        for _, row in sample_cov.iterrows():
            aa = row['amino_acid']
            if aa not in aa_index:
                continue
            pos = int(row['position']) - 1  # align_5p_idx is 1-indexed
            tlen = int(row['max_tRNA_len'])
            if pos < 0 or pos >= tlen:
                continue
            mapped_pos = len_map[tlen][pos]
            cov_count[aa_index[aa], mapped_pos] += row['count']

        # Cumulate 5'→3' (left to right) to get coverage
        for i in range(n_aa):
            for j in range(1, global_max_len):
                cov_count[i, j] += cov_count[i, j - 1]

        # Optional AA normalization: each AA equally weighted at 3' end
        if aa_norm:
            for i in range(n_aa):
                three_prime = cov_count[i, -1]
                if three_prime > 0:
                    cov_count[i, :] = cov_count[i, :] / three_prime / n_aa

        # Column-wise cumulation for stacked plots
        cov_count_sum = cov_count.copy()
        for i in range(1, n_aa):
            cov_count_sum[i] += cov_count_sum[i - 1]
        # Prepend zeros baseline row
        cov_count_sum = np.vstack((
            np.zeros(global_max_len), cov_count_sum
        ))

        return cov_count, cov_count_sum, aa_ordered, global_max_len

    # ---- panel 5: Behrens coverage plot ----

    def _panel_behrens_coverage(self):
        """Stacked step-area chart of read coverage by amino acid (Behrens-style).

        Matches the original TRNA_plot.plot_coverage(plot_type='behrens'):
        stacked filled step functions, each amino acid layered bottom-to-top.
        """
        if self.coverage_df is None or self.coverage_df.empty:
            return ''

        required = {'amino_acid', 'position', 'count', 'max_tRNA_len',
                     'sample_name_unique'}
        if not required.issubset(self.coverage_df.columns):
            return ''

        samples = sorted(self.coverage_df['sample_name_unique'].unique())
        if not samples:
            return ''

        aa_colors = _aa_palette()

        fig = go.Figure()
        buttons = []

        # Pre-compute matrices for all samples
        sample_data = []
        for sample in samples:
            scov = self.coverage_df[
                self.coverage_df['sample_name_unique'] == sample
            ]
            result = self._build_coverage_matrices(scov)
            sample_data.append(result)

        # Find consistent n_aas across samples for visibility toggling
        # Each sample produces: top outline + n_aa filled layers = n_aa + 1 traces
        traces_per_sample = []

        for s_idx, (sample, (cov_count, cov_count_sum, aa_ordered, max_len)) in enumerate(
            zip(samples, sample_data)
        ):
            if cov_count is None:
                traces_per_sample.append(0)
                continue

            visible = (s_idx == 0)
            n_aa = len(aa_ordered)
            x_vals = list(range(max_len))

            # Stacked step areas: draw each AA layer as fill between
            # cov_count_sum[i-1] (baseline) and cov_count_sum[i] (top)
            for i in range(n_aa):
                aa = aa_ordered[i]
                color = aa_colors.get(aa, '#b2bec3')
                top = cov_count_sum[i + 1]
                baseline = cov_count_sum[i]

                # Baseline trace (invisible, used as fill anchor)
                fig.add_trace(go.Scatter(
                    x=x_vals, y=baseline.tolist(),
                    mode='lines',
                    line=dict(width=0, color='rgba(0,0,0,0)'),
                    showlegend=False,
                    visible=visible,
                    hoverinfo='skip',
                    line_shape='hv',
                ))
                # Top trace with fill to previous
                fig.add_trace(go.Scatter(
                    x=x_vals, y=top.tolist(),
                    mode='lines',
                    name=aa,
                    line=dict(width=0.5, color=color),
                    line_shape='hv',
                    fill='tonexty',
                    fillcolor=color,
                    opacity=0.7,
                    visible=visible,
                    showlegend=(s_idx == 0),
                    hovertemplate=f'{aa}<br>'
                                 "Position: %{{x}}<br>"
                                 'Coverage: %{y:.0f}<extra></extra>',
                ))

            # Black outline on top
            fig.add_trace(go.Scatter(
                x=x_vals, y=cov_count_sum[-1].tolist(),
                mode='lines',
                line=dict(width=1, color='black'),
                line_shape='hv',
                showlegend=False,
                visible=visible,
                hoverinfo='skip',
            ))

            n_traces = n_aa * 2 + 1  # baseline+top per AA, plus outline
            traces_per_sample.append(n_traces)

            # Dropdown button
            total_traces = sum(traces_per_sample)
            vis = [False] * total_traces
            start = total_traces - n_traces
            for j in range(n_traces):
                vis[start + j] = True

        # Rebuild visibility arrays now that all traces are added
        total = sum(traces_per_sample)
        final_buttons = []
        offset = 0
        for s_idx, sample in enumerate(samples):
            n_t = traces_per_sample[s_idx]
            if n_t == 0:
                offset += n_t
                continue
            vis = [False] * total
            for j in range(n_t):
                vis[offset + j] = True
            final_buttons.append(dict(
                label=sample,
                method='update',
                args=[{'visible': vis}],
            ))
            offset += n_t

        if final_buttons:
            fig.update_layout(
                updatemenus=[dict(
                    type='dropdown',
                    direction='down',
                    x=1.0, xanchor='right',
                    y=1.15, yanchor='top',
                    buttons=final_buttons,
                    active=0,
                )],
            )

        fig.update_layout(
            title='tRNA Coverage by Amino Acid (Behrens Plot)',
            xaxis_title="5' to 3' index (mapped to longest tRNA)",
            yaxis_title='Read count',
            height=550,
        )
        div = self._fig_to_div(fig)
        return render_panel('Behrens Coverage Plot', div, anchor='frag-behrens')

    # ---- panel 6: needle coverage plot ----

    def _panel_needle_coverage(self):
        """Needle plot: symmetric funnel shapes per amino acid.

        Matches the original TRNA_plot.plot_coverage(plot_type='needle'):
        each AA is drawn as a symmetric area around its 3'-end midpoint,
        creating a distinctive "needle" shape that fans out from the center
        toward the 5' end proportional to coverage.
        """
        if self.coverage_df is None or self.coverage_df.empty:
            return ''

        required = {'amino_acid', 'position', 'count', 'max_tRNA_len',
                     'sample_name_unique'}
        if not required.issubset(self.coverage_df.columns):
            return ''

        samples = sorted(self.coverage_df['sample_name_unique'].unique())
        if not samples:
            return ''

        aa_colors = _aa_palette()

        fig = go.Figure()
        traces_per_sample = []

        # Pre-compute all sample data
        sample_data = []
        for sample in samples:
            scov = self.coverage_df[
                self.coverage_df['sample_name_unique'] == sample
            ]
            # Use aa_norm=True for needle plot (equal weight per AA at 3')
            result = self._build_coverage_matrices(scov, aa_norm=True)
            sample_data.append(result)

        for s_idx, (sample, (cov_count, cov_count_sum, aa_ordered, max_len)) in enumerate(
            zip(samples, sample_data)
        ):
            if cov_count is None:
                traces_per_sample.append(0)
                continue

            visible = (s_idx == 0)
            n_aa = len(aa_ordered)
            x_vals = list(range(max_len))

            # Compute 3'-end midpoints for each AA (from stacked cumulative)
            last_col = cov_count_sum[:, -1]  # shape (n_aa+1,)
            last_col_mid = np.zeros(n_aa)
            for i in range(n_aa):
                last_col_mid[i] = last_col[i] + (last_col[i + 1] - last_col[i]) / 2

            # Build top/bottom funnel curves for each AA
            cov_funnel_top = np.zeros((n_aa, max_len))
            cov_funnel_bot = np.zeros((n_aa, max_len))
            for j in range(max_len):
                cov = cov_count[:, j]
                cov_funnel_top[:, j] = last_col_mid + cov / 2
                cov_funnel_bot[:, j] = last_col_mid - cov / 2

            # Draw each needle as a filled area between bottom and top
            for i in range(n_aa):
                aa = aa_ordered[i]
                color = aa_colors.get(aa, '#b2bec3')

                # Bottom curve (invisible, fill anchor)
                fig.add_trace(go.Scatter(
                    x=x_vals, y=cov_funnel_bot[i].tolist(),
                    mode='lines',
                    line=dict(width=0, color='rgba(0,0,0,0)'),
                    showlegend=False,
                    visible=visible,
                    hoverinfo='skip',
                    line_shape='hv',
                ))
                # Top curve with fill to bottom
                fig.add_trace(go.Scatter(
                    x=x_vals, y=cov_funnel_top[i].tolist(),
                    mode='lines',
                    name=aa,
                    line=dict(width=0.5, color=color),
                    line_shape='hv',
                    fill='tonexty',
                    fillcolor=color,
                    opacity=0.7,
                    visible=visible,
                    showlegend=(s_idx == 0),
                    hovertemplate=f'{aa}<br>'
                                 "Position: %{{x}}<br>"
                                 'Coverage: %{y:.4f}<extra></extra>',
                ))

            n_traces = n_aa * 2  # baseline + top per AA
            traces_per_sample.append(n_traces)

        # Build dropdown buttons
        total = sum(traces_per_sample)
        final_buttons = []
        offset = 0
        for s_idx, sample in enumerate(samples):
            n_t = traces_per_sample[s_idx]
            if n_t == 0:
                offset += n_t
                continue
            vis = [False] * total
            for j in range(n_t):
                vis[offset + j] = True
            final_buttons.append(dict(
                label=sample,
                method='update',
                args=[{'visible': vis}],
            ))
            offset += n_t

        if final_buttons:
            fig.update_layout(
                updatemenus=[dict(
                    type='dropdown',
                    direction='down',
                    x=1.0, xanchor='right',
                    y=1.15, yanchor='top',
                    buttons=final_buttons,
                    active=0,
                )],
            )

        fig.update_layout(
            title='Per-AA Normalized Coverage (Needle Plot)',
            xaxis_title="5' to 3' index (mapped to longest tRNA)",
            yaxis_title="Normalized coverage\n(amino acids equally weighed at 3')",
            height=550,
        )
        div = self._fig_to_div(fig)
        return render_panel('Needle Coverage Plot', div, anchor='frag-needle')

    # ---- panel 7: synthetic tRNA integrity ----

    def _panel_synthetic_integrity(self):
        """Bar chart of integrity scores for synthetic control tRNAs."""
        if not self.source_prefixes:
            return ''

        df = self.frag_counts
        if df.empty or 'integrity_score' not in df.columns:
            return ''

        trna_col = 'tRNA_annotation' if 'tRNA_annotation' in df.columns else 'tRNA_name'
        if trna_col not in df.columns:
            return ''

        from trnaseq.charge.quantifier import classify_trna_source

        df = df.copy()
        df['_source'] = df[trna_col].apply(
            lambda x: classify_trna_source(x, self.source_prefixes)
        )
        syn = df[df['_source'] == 'synthetic']
        if syn.empty:
            return ''

        samples = sorted(syn['sample_name_unique'].unique())
        trnas = sorted(syn[trna_col].unique(), key=_isoacceptor_sort_key)
        palette = _discrete_palette(len(samples))

        fig = go.Figure()
        for i, sample in enumerate(samples):
            ssub = syn[syn['sample_name_unique'] == sample]
            vals = []
            for t in trnas:
                row = ssub[ssub[trna_col] == t]
                vals.append(float(row['integrity_score'].iloc[0])
                            if not row.empty else np.nan)
            short_names = ['-'.join(t.split('-')[1:]) if '-' in t else t
                           for t in trnas]
            fig.add_trace(go.Bar(
                x=short_names,
                y=vals,
                name=sample,
                marker_color=palette[i],
            ))

        fig.update_layout(
            barmode='group',
            title='Synthetic Control tRNA Integrity',
            xaxis_title='Synthetic tRNA',
            yaxis_title='Integrity score',
            yaxis_range=[0, 1.05],
            height=max(400, len(trnas) * 30),
        )
        div = self._fig_to_div(fig)
        return render_panel('Synthetic Control Integrity', div, anchor='frag-synthetic')
