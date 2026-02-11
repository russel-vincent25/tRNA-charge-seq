"""
tRNA-seq I/O Module

This module provides efficient data storage and loading utilities
using Parquet format for compressed, fast data access, plus batch
processing utilities for handling multiple samples.

Classes:
    tRNAseqDataStore: Main class for Parquet storage and retrieval
    BatchProcessor: Manager for parallel batch processing
    ProgressTracker: Progress tracking for batch operations

Functions:
    process_samples_parallel: Process multiple samples in parallel
    merge_sample_results: Merge results from multiple samples

Usage:
    from trnaseq.io import tRNAseqDataStore, BatchProcessor

    # Storage operations
    store = tRNAseqDataStore('project_dir/')
    store.save_charge_data(charge_df)
    data = store.load_charge_data(samples=['sample001'])

    # Batch processing
    processor = BatchProcessor(store, num_workers=8)
    results = processor.process_batch(
        samples=['s1', 's2', 's3'],
        func=process_sample_func
    )
"""

from .storage import tRNAseqDataStore
from .batch import BatchProcessor, ProgressTracker, process_samples_parallel, merge_sample_results

__all__ = [
    'tRNAseqDataStore',
    'BatchProcessor',
    'ProgressTracker',
    'process_samples_parallel',
    'merge_sample_results',
]
