"""
Fragment analysis module for tRNA-charge-seq.

Classifies reads into fragment types (full-length, RT drop-off, 5' tRF,
degraded), profiles RT drop-off positions, and computes fragment length
distributions from per-read stats CSVs.
"""

from trnaseq.fragments.analyser import FragmentAnalyser

__all__ = ['FragmentAnalyser']
