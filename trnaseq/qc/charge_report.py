"""
Charge QC Report Generator
============================

Generates charge_report.html — an interactive Plotly dashboard
summarizing tRNA charging levels across samples.

Panels:
    1. Per-sample mean charge bar chart
    2. Charge by amino acid heatmap
    3. Charge vs RPM scatter
    4. Replicate charge correlation
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional

import plotly.graph_objects as go
import plotly.offline
from plotly.subplots import make_subplots

from trnaseq.qc.report import _discrete_palette
from trnaseq.qc._common import (
    ReportContext,
    fig_to_div,
    host_filter,
    render_html_shell,
    render_panel,
)


class ChargeReportGenerator:
    """Generate interactive HTML dashboard for tRNA charge QC."""

    MAX_PAIRS = 30

    def __init__(
        self,
        charge_df_transcript: pd.DataFrame,
        charge_df_aa: pd.DataFrame,
        charge_summary: pd.DataFrame,
        sample_df: pd.DataFrame,
        source_prefixes: Optional[dict] = None,
        context: Optional[ReportContext] = None,
    ):
        """
        Parameters
        ----------
        charge_df_transcript : DataFrame
            Transcript-level charge data (from ChargeQuantifier).
        charge_df_aa : DataFrame
            Amino-acid-level charge data.
        charge_summary : DataFrame
            Summary statistics from ChargeQuantifier.
        sample_df : DataFrame
            Sample metadata with sample_name_unique, sample_name columns.
        source_prefixes : dict, optional
            {prefix: category} for classifying tRNA sources (e.g. synthetic).
        """
        self.charge_tr = charge_df_transcript.copy()
        self.charge_aa = charge_df_aa.copy()
        self.charge_summary = charge_summary.copy()
        self.sample_df = sample_df.copy()
        self.source_prefixes = source_prefixes
        self.context = context

    def generate_html_report(self, output_path) -> str:
        """Generate self-contained HTML report at *output_path*."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        panels = []

        p1 = self._panel_mean_charge_bar()
        if p1:
            panels.append(p1)

        p2 = self._panel_charge_aa_heatmap()
        if p2:
            panels.append(p2)

        p3 = self._panel_charge_vs_rpm()
        if p3:
            panels.append(p3)

        p4 = self._panel_replicate_charge_correlation()
        if p4:
            panels.append(p4)

        p5 = self._panel_synthetic_control_bar()
        if p5:
            panels.append(p5)

        body = '\n'.join(panels)

        html = render_html_shell(
            title='tRNA Charge Report',
            h1_text='tRNA Charge Dashboard',
            accent_color='#00b894',
            body=body,
            context=self.context,
        )

        output_path.write_text(html)
        return str(output_path)

    # ---- helper ----

    def _fig_to_div(self, fig):
        return fig_to_div(fig)

    def _host_filter(self, df):
        """Filter to host tRNAs if a source column is present.

        Thin wrapper around :func:`trnaseq.qc._common.host_filter` so existing
        call sites keep working after the helper was lifted module-level.
        """
        return host_filter(df)

    # ---- panel 1: mean charge bar ----

    def _panel_mean_charge_bar(self):
        """Per-sample mean charge bar chart."""
        df = self._host_filter(self.charge_tr)
        if df.empty or 'charge_canonical' not in df.columns:
            return ''

        per_sample = df.groupby('sample_name_unique')['charge_canonical'].agg(
            ['mean', 'std']).reset_index()
        per_sample.columns = ['sample_name_unique', 'mean_charge', 'std_charge']
        per_sample = per_sample.sort_values('mean_charge')

        colors = []
        for v in per_sample['mean_charge']:
            if v >= 50:
                colors.append('#00b894')
            elif v >= 30:
                colors.append('#fdcb6e')
            else:
                colors.append('#d63031')

        fig = go.Figure(go.Bar(
            x=per_sample['sample_name_unique'],
            y=per_sample['mean_charge'],
            error_y=dict(type='data', array=per_sample['std_charge'].fillna(0)),
            marker_color=colors,
            text=per_sample['mean_charge'].round(1),
            textposition='auto',
        ))
        fig.update_layout(
            title='Mean Charge per Sample (host tRNAs, canonical)',
            xaxis_title='Sample',
            yaxis_title='Mean charge (%)',
            yaxis_range=[0, 105],
            height=max(400, len(per_sample) * 25),
        )
        div = self._fig_to_div(fig)
        return render_panel('Mean Charge per Sample', div, anchor='charge-mean')

    # ---- panel 2: amino acid heatmap ----

    def _panel_charge_aa_heatmap(self):
        """Heatmap: amino acid x sample, color = charge_canonical."""
        df = self._host_filter(self.charge_aa)
        if df.empty or 'charge_canonical' not in df.columns:
            return ''

        pivot = df.pivot_table(
            index='amino_acid', columns='sample_name_unique',
            values='charge_canonical', aggfunc='mean')

        if pivot.empty:
            return ''

        # Sort amino acids alphabetically
        pivot = pivot.sort_index()

        fig = go.Figure(go.Heatmap(
            z=pivot.values,
            x=pivot.columns.tolist(),
            y=pivot.index.tolist(),
            colorscale='RdYlGn',
            zmin=0, zmax=100,
            colorbar=dict(title='Charge %'),
            hovertemplate='AA: %{y}<br>Sample: %{x}<br>Charge: %{z:.1f}%<extra></extra>',
        ))
        fig.update_layout(
            title='Charge by Amino Acid (host tRNAs)',
            height=max(400, len(pivot) * 22),
        )
        div = self._fig_to_div(fig)
        return render_panel('Charge by Amino Acid', div, anchor='charge-aa')

    # ---- panel 3: charge vs RPM ----

    def _panel_charge_vs_rpm(self):
        """Scatter: log10(RPM+1) vs charge_canonical, colored by amino acid."""
        df = self._host_filter(self.charge_tr)
        if df.empty or 'RPM' not in df.columns or 'charge_canonical' not in df.columns:
            return ''

        df = df.dropna(subset=['charge_canonical', 'RPM']).copy()
        if df.empty:
            return ''

        df['log10_RPM'] = np.log10(df['RPM'] + 1)

        # Cap legend at 15 amino acids + Other
        aa_col = 'amino_acid' if 'amino_acid' in df.columns else None
        if aa_col is None:
            return ''

        top_aas = df[aa_col].value_counts().head(15).index.tolist()
        df['aa_group'] = df[aa_col].where(df[aa_col].isin(top_aas), 'Other')

        groups = sorted(df['aa_group'].unique())
        palette = _discrete_palette(len(groups))
        color_map = {g: palette[i] for i, g in enumerate(groups)}

        fig = go.Figure()
        for group in groups:
            sub = df[df['aa_group'] == group]
            fig.add_trace(go.Scatter(
                x=sub['log10_RPM'],
                y=sub['charge_canonical'],
                mode='markers',
                name=group,
                marker=dict(size=5, color=color_map[group], opacity=0.6),
                hovertemplate='%{text}<br>log10(RPM+1): %{x:.2f}<br>Charge: %{y:.1f}%',
                text=sub.get('tRNA_anno_short', sub.get('tRNA_annotation', '')),
            ))

        fig.update_layout(
            title='Charge vs Abundance (host tRNAs)',
            xaxis_title='log10(RPM + 1)',
            yaxis_title='Charge canonical (%)',
            yaxis_range=[0, 105],
            height=550,
        )
        div = self._fig_to_div(fig)
        return render_panel('Charge vs Abundance', div, anchor='charge-vs-rpm')

    # ---- panel 4: replicate correlation ----

    def _panel_replicate_charge_correlation(self):
        """Pairwise scatter of charge_canonical between replicates."""
        df = self._host_filter(self.charge_tr)
        if df.empty or 'charge_canonical' not in df.columns:
            return ''
        if 'sample_name' not in self.sample_df.columns:
            return ''

        # Pivot: tRNA_annotation x sample_name_unique -> charge_canonical
        pivot = df.pivot_table(
            index='tRNA_annotation', columns='sample_name_unique',
            values='charge_canonical', aggfunc='mean')

        if pivot.shape[1] < 2:
            return ''

        # Group by sample_name to find replicate pairs
        sn_map = dict(zip(self.sample_df['sample_name_unique'],
                          self.sample_df['sample_name'].astype(str)))

        groups = {}
        for col in pivot.columns:
            sn = str(sn_map.get(col, col))
            groups.setdefault(sn, []).append(col)

        pairs = []
        for sn, members in groups.items():
            if len(members) >= 2:
                for i in range(len(members)):
                    for j in range(i + 1, len(members)):
                        pairs.append((members[i], members[j], sn))

        if not pairs:
            return ''

        if len(pairs) > self.MAX_PAIRS:
            return (f'<div style="padding:20px"><h3>Replicate Charge Correlation</h3>'
                    f'<p>Skipped: {len(pairs)} replicate pairs exceeds display '
                    f'limit ({self.MAX_PAIRS}).</p></div>')

        n_pairs = len(pairs)
        cols = min(3, n_pairs)
        rows = (n_pairs + cols - 1) // cols

        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=[f'{a} vs {b}' for a, b, _ in pairs],
        )

        for idx, (s1, s2, sn) in enumerate(pairs):
            r, c = divmod(idx, cols)

            # Get shared tRNAs
            mask = pivot[[s1, s2]].dropna()
            x = mask[s1].values
            y = mask[s2].values

            if len(x) < 2:
                continue

            pearson_r = np.corrcoef(x, y)[0, 1]

            fig.add_trace(go.Scatter(
                x=x, y=y, mode='markers',
                marker=dict(size=4, opacity=0.5),
                name=f'{sn} (R={pearson_r:.3f})',
                showlegend=True,
            ), row=r + 1, col=c + 1)

            # Diagonal
            max_val = max(x.max(), y.max(), 1)
            fig.add_trace(go.Scatter(
                x=[0, max_val], y=[0, max_val],
                mode='lines', line=dict(dash='dash', color='gray'),
                showlegend=False,
            ), row=r + 1, col=c + 1)

            fig.update_xaxes(title_text=f'Charge % ({s1})', row=r + 1, col=c + 1)
            fig.update_yaxes(title_text=f'Charge % ({s2})', row=r + 1, col=c + 1)

            ax_suffix = '' if idx == 0 else str(idx + 1)
            fig.add_annotation(
                x=0.05, y=0.95,
                xref=f'x{ax_suffix} domain', yref=f'y{ax_suffix} domain',
                text=f'R = {pearson_r:.3f}', showarrow=False,
                font=dict(size=12, color='red'),
            )

        fig.update_layout(
            title='Replicate Charge Correlation',
            height=350 * rows,
        )
        div = self._fig_to_div(fig)
        return render_panel('Replicate Charge Correlation', div, anchor='charge-replicate')

    # ---- panel 5: synthetic control bar ----

    def _panel_synthetic_control_bar(self):
        """Bar chart of charge_canonical per synthetic control tRNA."""
        if not self.source_prefixes:
            return ''

        df = self.charge_tr
        if df.empty or 'charge_canonical' not in df.columns:
            return ''

        trna_col = 'tRNA_annotation' if 'tRNA_annotation' in df.columns else None
        if trna_col is None:
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
        trnas = sorted(syn[trna_col].unique())
        palette = _discrete_palette(len(samples))

        fig = go.Figure()
        for i, sample in enumerate(samples):
            ssub = syn[syn['sample_name_unique'] == sample]
            vals = []
            for t in trnas:
                row = ssub[ssub[trna_col] == t]
                vals.append(float(row['charge_canonical'].iloc[0])
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
            title='Synthetic Control tRNA Charge (expected ~0%)',
            xaxis_title='Synthetic tRNA',
            yaxis_title='Charge canonical (%)',
            yaxis_range=[0, 105],
            height=max(400, len(trnas) * 30),
        )
        div = self._fig_to_div(fig)
        return render_panel('Synthetic Control Charge', div, anchor='charge-synthetic')
