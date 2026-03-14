"""
Differential tRNA Abundance Analysis (pyDESeq2 wrapper)
========================================================

Builds a count matrix from ALL_stats_aggregate.csv at user-specified
aggregation level (transcript, codon, aa), then runs DESeq2 via pyDESeq2.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations


class DifferentialAbundance:
    """
    Differential abundance analysis for tRNA-seq data.

    Parameters
    ----------
    stats_csv : str or Path
        Path to ALL_stats_aggregate.csv (from stage 2).
    sample_df : str, Path, or DataFrame
        Path to sample_df.xlsx or a DataFrame with sample_name_unique
        and sample_name columns.
    level : str
        Aggregation level: 'transcript', 'codon', or 'aa'.
    control_group : str or None
        Name of the control group (must match a sample_name value).
        Defaults to the first sample_name group alphabetically.
    """

    LEVEL_COLS = {
        'transcript': 'tRNA_annotation',
        'codon': 'codon',
        'aa': 'amino_acid',
    }

    def __init__(self, stats_csv, sample_df, level='aa', control_group=None):
        self.stats_csv = Path(stats_csv)
        self.level = level

        if level not in self.LEVEL_COLS:
            raise ValueError(f"level must be one of {list(self.LEVEL_COLS)}, got '{level}'")

        self.feature_col = self.LEVEL_COLS[level]

        # Load stats
        self.stats_df = pd.read_csv(self.stats_csv, keep_default_na=False)

        # Load sample metadata
        if isinstance(sample_df, (str, Path)):
            self.sample_meta = pd.read_excel(sample_df)
        else:
            self.sample_meta = sample_df.copy()

        # Ensure sample_name column exists
        if 'sample_name' not in self.sample_meta.columns:
            raise ValueError("sample_df must have a 'sample_name' column for group assignment")

        # Build condition map: sample_name_unique -> sample_name (group)
        self.condition_map = dict(zip(
            self.sample_meta['sample_name_unique'],
            self.sample_meta['sample_name'].astype(str),
        ))

        # Set control group
        groups = sorted(set(self.condition_map.values()))
        self._control_explicit = (control_group is not None)
        if control_group is not None:
            if control_group not in groups:
                raise ValueError(f"control_group '{control_group}' not in groups: {groups}")
            self.control_group = control_group
        else:
            self.control_group = groups[0]

        # Build count matrix
        self.count_matrix = self._build_count_matrix()
        self.results = None

    def _build_count_matrix(self):
        """Build samples x features count matrix."""
        df = self.stats_df
        # Filter to only samples present in metadata
        known_samples = set(self.condition_map)
        df = df[df['sample_name_unique'].isin(known_samples)]

        pivot = df.pivot_table(
            index='sample_name_unique',
            columns=self.feature_col,
            values='count',
            aggfunc='sum',
            fill_value=0,
        )
        return pivot.astype(int)

    def run_deseq2(self):
        """
        Run DESeq2 analysis via pyDESeq2.

        Returns
        -------
        DataFrame
            DESeq2 results with columns: baseMean, log2FoldChange, lfcSE,
            stat, pvalue, padj.
        """
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError:
            raise ImportError(
                "pyDESeq2 is required for differential abundance analysis.\n"
                "Install with: pip install pydeseq2"
            )

        # Build metadata DataFrame for pyDESeq2
        meta = pd.DataFrame({
            'condition': [self.condition_map[s] for s in self.count_matrix.index]
        }, index=self.count_matrix.index)

        # Initialize and run DESeq2
        dds = DeseqDataSet(
            counts=self.count_matrix,
            metadata=meta,
            design_factors='condition',
        )
        dds.deseq2()

        # Extract results: each group vs control (explicit) or all pairwise
        groups = sorted(set(meta['condition']))

        all_results = []
        if self._control_explicit:
            # Each non-control group vs the designated control
            other_groups = [g for g in groups if g != self.control_group]
            pairs = [(g, self.control_group) for g in other_groups]
        else:
            # All pairwise combinations, both directions
            pairs = []
            for g1, g2 in combinations(groups, 2):
                pairs.append((g1, g2))

        for g1, g2 in pairs:
            stat_res = DeseqStats(
                dds,
                contrast=['condition', g1, g2],
            )
            stat_res.summary()
            res_df = stat_res.results_df.copy()
            res_df['comparison'] = f'{g1}_vs_{g2}'
            res_df['feature'] = res_df.index
            all_results.append(res_df)

        if all_results:
            self.results = pd.concat(all_results, ignore_index=True)
        else:
            self.results = pd.DataFrame()

        return self.results

    def volcano_plot(self, output_path=None, padj_threshold=0.05, lfc_threshold=1.0):
        """
        Generate an interactive Plotly volcano plot.

        Parameters
        ----------
        output_path : str or Path, optional
            If provided, saves the plot as an HTML file.
        padj_threshold : float
            Adjusted p-value threshold for significance.
        lfc_threshold : float
            Log2 fold change threshold for significance.

        Returns
        -------
        plotly.graph_objects.Figure
        """
        if self.results is None or self.results.empty:
            raise ValueError("Run run_deseq2() first")

        import plotly.graph_objects as go

        df = self.results.copy()
        df['-log10(padj)'] = -np.log10(df['padj'].clip(lower=1e-300))

        # Classify points
        sig_up = (df['padj'] < padj_threshold) & (df['log2FoldChange'] > lfc_threshold)
        sig_down = (df['padj'] < padj_threshold) & (df['log2FoldChange'] < -lfc_threshold)
        ns = ~(sig_up | sig_down)

        fig = go.Figure()
        for mask, name, color in [
            (ns, 'Not significant', '#b2bec3'),
            (sig_up, 'Up', '#d63031'),
            (sig_down, 'Down', '#0984e3'),
        ]:
            sub = df[mask]
            if sub.empty:
                continue
            fig.add_trace(go.Scatter(
                x=sub['log2FoldChange'],
                y=sub['-log10(padj)'],
                mode='markers',
                name=name,
                marker=dict(color=color, size=6, opacity=0.7),
                text=sub['feature'],
                hovertemplate='%{text}<br>LFC: %{x:.2f}<br>-log10(padj): %{y:.2f}',
            ))

        fig.add_hline(y=-np.log10(padj_threshold), line_dash='dash', line_color='gray')
        fig.add_vline(x=lfc_threshold, line_dash='dash', line_color='gray')
        fig.add_vline(x=-lfc_threshold, line_dash='dash', line_color='gray')

        fig.update_layout(
            title='Volcano Plot: Differential tRNA Abundance',
            xaxis_title='log2(Fold Change)',
            yaxis_title='-log10(adjusted p-value)',
            height=550,
        )

        if output_path:
            import plotly.offline
            plotly.offline.plot(fig, filename=str(output_path), auto_open=False)

        return fig

    def ma_plot(self, output_path=None, padj_threshold=0.05):
        """
        Generate an interactive Plotly MA plot.

        Parameters
        ----------
        output_path : str or Path, optional
            If provided, saves the plot as an HTML file.
        padj_threshold : float
            Adjusted p-value threshold for significance.

        Returns
        -------
        plotly.graph_objects.Figure
        """
        if self.results is None or self.results.empty:
            raise ValueError("Run run_deseq2() first")

        import plotly.graph_objects as go

        df = self.results.copy()
        sig = df['padj'] < padj_threshold

        fig = go.Figure()
        for mask, name, color in [
            (~sig, 'Not significant', '#b2bec3'),
            (sig, 'Significant', '#d63031'),
        ]:
            sub = df[mask]
            if sub.empty:
                continue
            fig.add_trace(go.Scatter(
                x=sub['baseMean'],
                y=sub['log2FoldChange'],
                mode='markers',
                name=name,
                marker=dict(color=color, size=6, opacity=0.7),
                text=sub['feature'],
                hovertemplate='%{text}<br>baseMean: %{x:.1f}<br>LFC: %{y:.2f}',
            ))

        fig.add_hline(y=0, line_color='black', opacity=0.3)
        fig.update_layout(
            title='MA Plot: Differential tRNA Abundance',
            xaxis_title='Mean of normalized counts',
            xaxis_type='log',
            yaxis_title='log2(Fold Change)',
            height=550,
        )

        if output_path:
            import plotly.offline
            plotly.offline.plot(fig, filename=str(output_path), auto_open=False)

        return fig

    def export_results(self, output_dir):
        """
        Export all results to output directory.

        Saves:
            - count_matrix.csv
            - deseq2_results.csv
            - volcano.html
            - ma_plot.html
        """
        if self.results is None:
            raise ValueError("Run run_deseq2() first")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.count_matrix.to_csv(output_dir / 'count_matrix.csv')
        self.results.to_csv(output_dir / 'deseq2_results.csv', index=False)
        self.volcano_plot(output_dir / 'volcano.html')
        self.ma_plot(output_dir / 'ma_plot.html')
