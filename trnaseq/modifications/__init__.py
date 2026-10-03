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
- load_channel_priors: Per-RT-enzyme detection-channel priors
"""

from .rt_signatures import RTSignatureAnalyzer, analyze_rt_signatures
from .modification_caller import (
    ModificationCaller,
    ModificationProfile,
    MODIFICATION_PROFILES,
    benjamini_hochberg_fdr,
    ChannelBackground,
    estimate_background_error_rate,
    estimate_channel_backgrounds,
    ReplicateAggregator,
)
from .channel_priors import load_channel_priors, normalize_rt_enzyme
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
    'ChannelBackground',
    'estimate_channel_backgrounds',
    'ReplicateAggregator',
    'load_channel_priors',
    'normalize_rt_enzyme',
    'PositionalExtractor',
    'MODOMICSAnnotator',
]
