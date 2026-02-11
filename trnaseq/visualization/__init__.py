"""
tRNA-seq Visualization Module

This module provides visualization tools for tRNA-seq data,
including alignment viewers and coverage plots.

Classes:
    AlignmentViewer: Interactive alignment viewer for JSON.bz2 files

Usage:
    from trnaseq.visualization import AlignmentViewer

    viewer = AlignmentViewer('sample001_SWalign.json.bz2')
    viewer.quick_view('tRNA-Ala-TGC-1')
"""

from .alignment_viewer import AlignmentViewer

__all__ = ['AlignmentViewer']
