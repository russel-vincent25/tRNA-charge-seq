"""
QC Summary Report Generator
============================

Generates QC_summary.csv and QC_report.html from pipeline outputs.
No new computation — reads existing DataFrames produced by earlier stages.

Outputs:
    QC_summary.csv  — one row per sample, all pipeline metrics
    QC_report.html  — interactive Plotly dashboard (self-contained HTML)
"""

import pandas as pd
import numpy as np
from pathlib import Path

import plotly.graph_objects as go
import plotly.offline
from plotly.subplots import make_subplots


class QCReportGenerator:
    """Generate QC summary CSV and HTML dashboard from pipeline outputs."""

    # Threshold definitions: (warn_below, fail_below) — higher is better
    THRESHOLDS = {
        'percent_successfully_merged': (90, 70),
        'percent_BC-mapped':           (75, 50),
        'percent_seqs_after_UMI_trim': (90, 70),
        'Mapping_percent':             (70, 40),
        'percent_single_annotation':   (80, 60),
    }

    def __init__(self, project_dir, sample_df, inp_file_df,
                 charge_summary_df=None, stats_df=None, bc_dir=None):
        """
        Parameters
        ----------
        project_dir : Path
            Project root directory.
        sample_df : DataFrame
            Sample info with QC metrics (from stages 0–2).
        inp_file_df : DataFrame
            Input-file-level metrics (merge rate, BC mapping).
        charge_summary_df : DataFrame, optional
            Charge summary (from stage 3). None if stage 3 was skipped.
        stats_df : DataFrame, optional
            ALL_stats_aggregate for composition analysis. None if unavailable.
        bc_dir : Path, optional
            BC_split directory containing read_length_distributions.csv.
        """
        self.project_dir = Path(project_dir)
        self.sample_df = sample_df.copy()
        self.inp_file_df = inp_file_df.copy()
        self.charge_summary_df = charge_summary_df
        self.stats_df = stats_df
        self.bc_dir = Path(bc_dir) if bc_dir is not None else None

    # ------------------------------------------------------------------
    # QC_summary.csv
    # ------------------------------------------------------------------

    def build_summary(self):
        """Build QC summary DataFrame — one row per sample."""

        df = self.sample_df.copy()

        # Merge inp_file_df metrics (merge %, BC mapping %) onto sample rows
        if not self.inp_file_df.empty:
            inp_cols = ['fastq_mate1_filename',
                        'N_pairs', 'N_merged', 'percent_successfully_merged',
                        'percent_BC-mapped']
            available = [c for c in inp_cols if c in self.inp_file_df.columns]
            if available and 'fastq_mate1_filename' in available:
                df = df.merge(
                    self.inp_file_df[available].drop_duplicates(),
                    on='fastq_mate1_filename', how='left', suffixes=('', '_inp')
                )

        # Add charge mean per sample (if stage 3 ran)
        if self.charge_summary_df is not None and len(self.charge_summary_df) > 0:
            charge_file = self.project_dir / 'charge_analysis' / 'charge_df_transcript.csv'
            if charge_file.exists():
                charge_df = pd.read_csv(charge_file)
                if 'charge_canonical' in charge_df.columns:
                    per_sample = (charge_df
                                  .groupby('sample_name_unique')['charge_canonical']
                                  .mean()
                                  .reset_index()
                                  .rename(columns={'charge_canonical': 'charge_canonical_mean'}))
                    df = df.merge(per_sample, on='sample_name_unique', how='left')

        # Select columns for QC summary (keep what exists)
        desired = [
            'sample_name_unique', 'sample_name', 'species',
            'N_pairs', 'N_merged', 'percent_successfully_merged',
            'percent_BC-mapped',
            'N_total', 'N_after_trim',
            'percent_seqs_after_UMI_trim', 'percent_UMI_obs-vs-exp',
            'N_UMI_observed', 'N_UMI_expected',
            'Mapping_percent', 'percent_single_annotation',
            'N_full_length', 'N_rt_dropoff', 'N_5p_fragment', 'N_degraded',
            'N_total_aligned',
            'charge_canonical_mean',
        ]
        out_cols = [c for c in desired if c in df.columns]
        return df[out_cols]

    def save_summary_csv(self, output_path=None):
        """Write QC_summary.csv and return the DataFrame."""
        summary = self.build_summary()
        if output_path is None:
            qc_dir = self.project_dir / 'qc_reports'
            qc_dir.mkdir(exist_ok=True)
            output_path = qc_dir / 'QC_summary.csv'
        summary.to_csv(output_path, index=False)
        return summary

    # ------------------------------------------------------------------
    # QC_report.html
    # ------------------------------------------------------------------

    def generate_html_report(self, output_path=None):
        """Generate self-contained interactive HTML dashboard."""
        summary = self.build_summary()
        if output_path is None:
            qc_dir = self.project_dir / 'qc_reports'
            qc_dir.mkdir(exist_ok=True)
            output_path = qc_dir / 'QC_report.html'

        panels = []

        # 1. Summary table
        panels.append(self._html_summary_table(summary))

        # 2. Pipeline funnel
        funnel = self._plotly_funnel(summary)
        if funnel:
            panels.append(funnel)

        # 3. Mapping rate
        if 'Mapping_percent' in summary.columns:
            panels.append(self._plotly_mapping_bar(summary))

        # 4. UMI saturation
        if all(c in summary.columns for c in ['N_UMI_observed', 'N_UMI_expected']):
            panels.append(self._plotly_umi_saturation(summary))

        # 5. Fragment types
        frag_cols = ['N_full_length', 'N_rt_dropoff', 'N_5p_fragment', 'N_degraded']
        if any(c in summary.columns for c in frag_cols):
            panels.append(self._plotly_fragment_bar(summary))

        # 6. Read length distribution
        readlen = self._plotly_read_length_density()
        if readlen:
            panels.append(readlen)

        # 7. PCA
        if self.stats_df is not None:
            pca = self._plotly_pca()
            if pca:
                panels.append(pca)

        # 8. Replicate correlation
        if self.stats_df is not None:
            corr = self._plotly_replicate_correlation()
            if corr:
                panels.append(corr)

        body = '\n'.join(panels)

        # Embed Plotly.js once
        plotlyjs = plotly.offline.get_plotlyjs()

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>tRNA-charge-seq QC Report</title>
<script>{plotlyjs}</script>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
       margin: 0; padding: 20px 40px; background: #f5f6fa; color: #2d3436; }}
h1 {{ border-bottom: 3px solid #0984e3; padding-bottom: 10px; }}
h2 {{ color: #2d3436; margin-top: 40px; }}
.card {{ background: white; border-radius: 8px; padding: 20px; margin: 20px 0;
         box-shadow: 0 2px 8px rgba(0,0,0,0.08); overflow-x: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
th {{ background: #0984e3; color: white; padding: 8px 12px; text-align: left;
      position: sticky; top: 0; }}
td {{ padding: 6px 12px; border-bottom: 1px solid #dfe6e9; }}
tr:hover {{ background: #f0f3f8; }}
.pass {{ background-color: #00b89433; }}
.warn {{ background-color: #fdcb6e55; }}
.fail {{ background-color: #d6336c33; }}
.footer {{ font-size: 11px; color: #636e72; margin-top: 40px; text-align: center; }}
</style>
</head>
<body>
<h1>tRNA-charge-seq QC Dashboard</h1>
<p>Project: <code>{self.project_dir}</code></p>
{body}
<div class="footer">Generated by trnaseq QC module</div>
</body>
</html>"""

        Path(output_path).write_text(html)
        return output_path

    # ---- helper: embed a plotly figure as div ----

    def _fig_to_div(self, fig):
        """Convert a Plotly figure to an HTML div string (no plotlyjs)."""
        return plotly.offline.plot(fig, output_type='div', include_plotlyjs=False)

    # ---- panel: summary table ----

    def _html_summary_table(self, summary):
        """Render summary table with pass/warn/fail coloring."""
        header_cells = ''.join(f'<th>{c}</th>' for c in summary.columns)
        rows = []
        for _, row in summary.iterrows():
            cells = []
            for col in summary.columns:
                val = row[col]
                css = ''
                if col in self.THRESHOLDS and pd.notna(val):
                    warn, fail = self.THRESHOLDS[col]
                    if val < fail:
                        css = ' class="fail"'
                    elif val < warn:
                        css = ' class="warn"'
                    else:
                        css = ' class="pass"'
                fmt_val = f'{val:.1f}' if isinstance(val, float) else str(val)
                cells.append(f'<td{css}>{fmt_val}</td>')
            rows.append('<tr>' + ''.join(cells) + '</tr>')

        return f"""
<h2>Sample Summary</h2>
<div class="card">
<table>
<tr>{header_cells}</tr>
{''.join(rows)}
</table>
</div>"""

    # ---- panel: pipeline funnel ----

    def _plotly_funnel(self, summary):
        """Pipeline funnel: reads surviving each stage."""
        count_cols = [c for c in ['N_pairs', 'N_merged', 'N_total', 'N_after_trim']
                      if c in summary.columns]
        if len(count_cols) < 2:
            return ''

        labels = {
            'N_pairs': 'Raw pairs',
            'N_merged': 'Merged',
            'N_total': 'After BC split',
            'N_after_trim': 'After UMI trim',
        }

        fig = go.Figure()
        for _, row in summary.iterrows():
            vals = [row.get(c, 0) for c in count_cols]
            fig.add_trace(go.Scatter(
                x=[labels.get(c, c) for c in count_cols],
                y=vals,
                mode='lines+markers',
                name=row.get('sample_name_unique', ''),
            ))
        fig.update_layout(
            title='Pipeline Funnel: Reads Surviving Each Stage',
            yaxis_title='Read count',
            hovermode='x unified',
            height=450,
        )
        div = self._fig_to_div(fig)
        return f'<h2>Pipeline Funnel</h2><div class="card">{div}</div>'

    # ---- panel: mapping rate bar ----

    def _plotly_mapping_bar(self, summary):
        """Per-sample mapping rate horizontal bar chart."""
        ordered = summary.sort_values('Mapping_percent')
        colors = ['#d63031' if v < 40 else '#fdcb6e' if v < 70 else '#00b894'
                  for v in ordered['Mapping_percent']]

        fig = go.Figure(go.Bar(
            x=ordered['Mapping_percent'],
            y=ordered['sample_name_unique'],
            orientation='h',
            marker_color=colors,
            text=ordered['Mapping_percent'].round(1),
            textposition='auto',
        ))
        fig.add_vline(x=70, line_dash='dash', line_color='gray',
                      annotation_text='warn=70%')
        fig.update_layout(
            title='SWIPE Mapping Rate per Sample',
            xaxis_title='Mapping %',
            xaxis_range=[0, 105],
            height=max(300, len(summary) * 40),
        )
        div = self._fig_to_div(fig)
        return f'<h2>Mapping Rate</h2><div class="card">{div}</div>'

    # ---- panel: UMI saturation ----

    def _plotly_umi_saturation(self, summary):
        """UMI observed vs expected scatter."""
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=summary['N_UMI_expected'],
            y=summary['N_UMI_observed'],
            mode='markers+text',
            text=summary['sample_name_unique'],
            textposition='top right',
            textfont_size=9,
            marker=dict(size=10),
            name='Samples',
        ))
        max_val = max(summary['N_UMI_expected'].max(), summary['N_UMI_observed'].max()) * 1.1
        fig.add_trace(go.Scatter(
            x=[0, max_val], y=[0, max_val],
            mode='lines', line=dict(dash='dash', color='black'),
            opacity=0.3, name='y=x (full complexity)',
        ))
        fig.update_layout(
            title='UMI Saturation (Library Complexity)',
            xaxis_title='Expected UMIs',
            yaxis_title='Observed UMIs',
            height=500,
        )
        div = self._fig_to_div(fig)
        return f'<h2>UMI Saturation</h2><div class="card">{div}</div>'

    # ---- panel: fragment types ----

    def _plotly_fragment_bar(self, summary):
        """Stacked bar of fragment types per sample."""
        frag_cols = ['N_full_length', 'N_rt_dropoff', 'N_5p_fragment', 'N_degraded']
        present = [c for c in frag_cols if c in summary.columns]
        if not present:
            return ''

        colors = {
            'N_full_length': '#00b894',
            'N_rt_dropoff': '#0984e3',
            'N_5p_fragment': '#fdcb6e',
            'N_degraded': '#d63031',
        }
        labels = {
            'N_full_length': 'Full-length',
            'N_rt_dropoff': 'RT drop-off',
            'N_5p_fragment': "5' fragment",
            'N_degraded': 'Degraded',
        }

        fig = go.Figure()
        for col in present:
            fig.add_trace(go.Bar(
                x=summary['sample_name_unique'],
                y=summary[col],
                name=labels.get(col, col),
                marker_color=colors.get(col, '#b2bec3'),
            ))
        fig.update_layout(
            barmode='stack',
            title='Fragment Type Distribution per Sample',
            yaxis_title='Read count',
            height=450,
        )
        div = self._fig_to_div(fig)
        return f'<h2>Fragment Types</h2><div class="card">{div}</div>'

    # ---- panel: read length density ----

    def _plotly_read_length_density(self):
        """Overlapped line density plot of read lengths from BC_split."""
        if self.bc_dir is None:
            return ''
        csv_path = self.bc_dir / 'read_length_distributions.csv'
        if not csv_path.exists():
            return ''

        rl_df = pd.read_csv(csv_path)
        if rl_df.empty:
            return ''

        fig = go.Figure()
        for sample in sorted(rl_df['sample_name_unique'].unique()):
            sub = rl_df[rl_df['sample_name_unique'] == sample].sort_values('read_length')
            fig.add_trace(go.Scatter(
                x=sub['read_length'],
                y=sub['count'],
                mode='lines',
                name=sample,
            ))
        fig.update_layout(
            title='Read Length Distribution (Post-BC Split)',
            xaxis_title='Read length (nt)',
            yaxis_title='Count',
            hovermode='x unified',
            height=450,
        )
        div = self._fig_to_div(fig)
        return f'<h2>Read Length Distribution</h2><div class="card">{div}</div>'

    # ---- panel: PCA ----

    def _plotly_pca(self):
        """PCA of samples based on tRNA RPM profiles (numpy SVD)."""
        df = self.stats_df
        if df is None or 'count' not in df.columns:
            return ''

        # Build RPM matrix: samples x tRNA annotations
        sample_col = 'sample_name_unique'
        trna_col = 'tRNA_annotation'
        if sample_col not in df.columns or trna_col not in df.columns:
            return ''

        pivot = df.pivot_table(
            index=sample_col, columns=trna_col,
            values='count', aggfunc='sum', fill_value=0)

        if pivot.shape[0] < 3:
            return ''

        # RPM normalization
        totals = pivot.sum(axis=1)
        rpm = pivot.div(totals, axis=0) * 1e6

        # Center and SVD
        X = rpm.values
        X_centered = X - X.mean(axis=0)
        try:
            U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
        except np.linalg.LinAlgError:
            return ''

        pc1 = U[:, 0] * S[0]
        pc2 = U[:, 1] * S[1]
        var_explained = (S ** 2) / (S ** 2).sum() * 100

        # Color by sample_name (replicate group) if available
        samples = pivot.index.tolist()
        if 'sample_name' in self.sample_df.columns:
            sn_map = dict(zip(self.sample_df['sample_name_unique'], self.sample_df['sample_name']))
            groups = [str(sn_map.get(s, s)) for s in samples]
        else:
            groups = samples

        unique_groups = list(dict.fromkeys(groups))
        palette = _discrete_palette(len(unique_groups))
        color_map = {g: palette[i] for i, g in enumerate(unique_groups)}

        fig = go.Figure()
        for group in unique_groups:
            mask = [i for i, g in enumerate(groups) if g == group]
            fig.add_trace(go.Scatter(
                x=[pc1[i] for i in mask],
                y=[pc2[i] for i in mask],
                mode='markers+text',
                text=[samples[i] for i in mask],
                textposition='top right',
                textfont_size=9,
                marker=dict(size=10, color=color_map[group]),
                name=group,
            ))
        fig.update_layout(
            title='PCA of tRNA RPM Profiles',
            xaxis_title=f'PC1 ({var_explained[0]:.1f}%)',
            yaxis_title=f'PC2 ({var_explained[1]:.1f}%)',
            height=500,
        )
        div = self._fig_to_div(fig)
        return f'<h2>PCA</h2><div class="card">{div}</div>'

    # ---- panel: replicate correlation ----

    def _plotly_replicate_correlation(self):
        """Log-scale scatter + Pearson R for replicate pairs."""
        df = self.stats_df
        if df is None or 'sample_name' not in self.sample_df.columns:
            return ''

        sample_col = 'sample_name_unique'
        trna_col = 'tRNA_annotation'
        if sample_col not in df.columns or trna_col not in df.columns:
            return ''

        # Build RPM matrix
        pivot = df.pivot_table(
            index=sample_col, columns=trna_col,
            values='count', aggfunc='sum', fill_value=0)

        totals = pivot.sum(axis=1)
        rpm = pivot.div(totals, axis=0) * 1e6

        # Group by sample_name to find replicate pairs
        sn_map = dict(zip(self.sample_df['sample_name_unique'], self.sample_df['sample_name']))
        groups = {}
        for snu in rpm.index:
            sn = str(sn_map.get(snu, snu))
            groups.setdefault(sn, []).append(snu)

        pairs = []
        for sn, members in groups.items():
            if len(members) >= 2:
                for i in range(len(members)):
                    for j in range(i + 1, len(members)):
                        pairs.append((members[i], members[j], sn))

        if not pairs:
            return ''

        MAX_PAIRS = 30
        if len(pairs) > MAX_PAIRS:
            return (f'<div style="padding:20px"><h3>Replicate Correlation</h3>'
                    f'<p>Skipped: {len(pairs)} replicate pairs exceeds display '
                    f'limit ({MAX_PAIRS}). Use PCA panel for replicate QC.</p></div>')

        n_pairs = len(pairs)
        cols = min(3, n_pairs)
        rows = (n_pairs + cols - 1) // cols

        fig = make_subplots(
            rows=rows, cols=cols,
            subplot_titles=[f'{a} vs {b}' for a, b, _ in pairs],
        )

        for idx, (s1, s2, sn) in enumerate(pairs):
            r, c = divmod(idx, cols)
            x = np.log10(rpm.loc[s1].values + 1)
            y = np.log10(rpm.loc[s2].values + 1)

            pearson_r = np.corrcoef(x, y)[0, 1]

            fig.add_trace(go.Scatter(
                x=x, y=y, mode='markers',
                marker=dict(size=3, opacity=0.5),
                name=f'{sn} (R={pearson_r:.3f})',
                showlegend=True,
            ), row=r + 1, col=c + 1)

            # Add diagonal
            max_val = max(x.max(), y.max())
            fig.add_trace(go.Scatter(
                x=[0, max_val], y=[0, max_val],
                mode='lines', line=dict(dash='dash', color='gray'),
                showlegend=False,
            ), row=r + 1, col=c + 1)

            fig.update_xaxes(title_text=f'log10(RPM+1) {s1}', row=r + 1, col=c + 1)
            fig.update_yaxes(title_text=f'log10(RPM+1) {s2}', row=r + 1, col=c + 1)

            # Annotate Pearson R
            # Plotly uses 'x domain' for first subplot, 'x2 domain' for second, etc.
            ax_suffix = '' if idx == 0 else str(idx + 1)
            fig.add_annotation(
                x=0.05, y=0.95,
                xref=f'x{ax_suffix} domain', yref=f'y{ax_suffix} domain',
                text=f'R = {pearson_r:.3f}', showarrow=False,
                font=dict(size=12, color='red'),
            )

        fig.update_layout(
            title='Replicate Correlation (log10 RPM)',
            height=350 * rows,
        )
        div = self._fig_to_div(fig)
        return f'<h2>Replicate Correlation</h2><div class="card">{div}</div>'


def _discrete_palette(n):
    """Return n distinct hex colors from matplotlib tab20."""
    import matplotlib
    cmap = matplotlib.colormaps['tab20']
    tab20 = [matplotlib.colors.rgb2hex(cmap(i)) for i in range(20)]
    if n <= 20:
        return tab20[:n]
    # Extend with evenly spaced hues for n > 20
    import colorsys
    extra = []
    for i in range(n - 20):
        h = (20 + i) / (n + 1)
        r, g, b = colorsys.hls_to_rgb(h, 0.5, 0.7)
        extra.append(f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}')
    return tab20 + extra
