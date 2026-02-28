"""
tRNA Modification Analysis Module

This module provides tools for detecting and analyzing tRNA modifications
based on RT signature analysis from sequencing data.

Components:
- PositionalExtractor: Stream SWalign JSONs to build per-position count matrices
- RTSignatureAnalyzer: Analyze mismatch, gap, and RT stop signatures
- ModificationCaller: Call known and novel modifications from signatures
- MODOMICSAnnotator: Integrate MODOMICS database for modification annotation
- ModificationProfile: Dataclass defining known modification RT signatures
"""

from .rt_signatures import RTSignatureAnalyzer, analyze_rt_signatures
from .modification_caller import (
    ModificationCaller,
    ModificationProfile,
    MODIFICATION_PROFILES,
    benjamini_hochberg_fdr,
    estimate_background_error_rate,
    ReplicateAggregator,
)
from .positional import PositionalExtractor
from .modomics import MODOMICSAnnotator

__all__ = [
    'RTSignatureAnalyzer',
    'analyze_rt_signatures',
    'ModificationCaller',
    'ModificationProfile',
    'MODIFICATION_PROFILES',
    'benjamini_hochberg_fdr',
    'estimate_background_error_rate',
    'ReplicateAggregator',
    'PositionalExtractor',
    'MODOMICSAnnotator',
]
