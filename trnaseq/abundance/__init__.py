"""
tRNA Abundance & Differential Expression Module
================================================

Wraps pyDESeq2 for differential tRNA abundance analysis.

Usage:
    from trnaseq.abundance import DifferentialAbundance

    da = DifferentialAbundance(
        stats_csv='ALL_stats_aggregate.csv',
        sample_df='sample_df.xlsx',
        level='aa',
        control_group='WT',
    )
    results = da.run_deseq2()
    da.volcano_plot('results/volcano.html')
"""

from trnaseq.abundance.differential import DifferentialAbundance

__all__ = ['DifferentialAbundance']
