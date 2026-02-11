"""
tRNA Charge Quantification Module

This module provides tools for quantifying tRNA charging levels
from tRNA-seq alignment statistics.

Classes:
    ChargeQuantifier: Main class for charge quantification from CSV stats

Usage:
    from trnaseq.charge import ChargeQuantifier

    quantifier = ChargeQuantifier(stats_csv='ALL_stats_aggregate.csv')
    charge_df = quantifier.quantify_all()
"""

from .quantifier import ChargeQuantifier

__all__ = ['ChargeQuantifier']
