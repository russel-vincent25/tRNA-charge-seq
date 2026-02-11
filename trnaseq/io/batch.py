"""
Batch processing utilities for efficient multi-sample data processing.

This module provides utilities for processing multiple samples in parallel,
tracking progress, and managing batch operations on tRNA-seq data.

Classes:
    BatchProcessor: Manager for parallel batch processing
    ProgressTracker: Progress tracking for batch operations

Functions:
    process_samples_parallel: Process multiple samples in parallel
    merge_sample_results: Merge results from multiple samples
"""

import os
import time
from pathlib import Path
from typing import Callable, List, Optional, Dict, Any, Union
from multiprocessing import Pool
import logging

import pandas as pd
from .storage import tRNAseqDataStore


class ProgressTracker:
    """
    Simple progress tracker for batch operations.

    Tracks completion count, timing, and provides formatted progress output.

    Example:
        >>> tracker = ProgressTracker(total=256)
        >>> for i in range(256):
        ...     # process item
        ...     tracker.update(f"Processing sample {i}")
    """

    def __init__(self, total: int, desc: str = "Processing"):
        """
        Initialize progress tracker.

        Args:
            total: Total number of items to process
            desc: Description for logging
        """
        self.total = total
        self.desc = desc
        self.completed = 0
        self.start_time = time.time()
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """Setup logger for progress tracking."""
        logger = logging.getLogger(f"ProgressTracker-{self.desc}")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger

    def update(self, message: str = "") -> None:
        """
        Update progress.

        Args:
            message: Optional message to log
        """
        self.completed += 1
        elapsed = time.time() - self.start_time
        rate = self.completed / elapsed if elapsed > 0 else 0
        remaining = (self.total - self.completed) / rate if rate > 0 else 0

        progress_str = (
            f"[{self.desc}] {self.completed}/{self.total} "
            f"({100*self.completed/self.total:.1f}%) "
            f"| {rate:.2f} items/sec | ETA: {remaining:.0f}s"
        )

        if message:
            progress_str += f" | {message}"

        self.logger.info(progress_str)

    def summary(self) -> None:
        """Print completion summary."""
        elapsed = time.time() - self.start_time
        rate = self.completed / elapsed if elapsed > 0 else 0

        summary_str = (
            f"[{self.desc}] COMPLETED: {self.completed}/{self.total} "
            f"in {elapsed:.1f}s ({rate:.2f} items/sec)"
        )
        self.logger.info(summary_str)


class BatchProcessor:
    """
    Manager for parallel batch processing of samples.

    Handles distributed processing of multiple samples with progress tracking,
    result aggregation, and error handling.

    Attributes:
        store: tRNAseqDataStore instance
        num_workers: Number of parallel workers
        chunk_size: Batch size for processing

    Example:
        >>> processor = BatchProcessor(store, num_workers=8)
        >>> results = processor.process_batch(
        ...     samples=['s1', 's2', 's3'],
        ...     func=process_sample,
        ...     output_type='charge'
        ... )
    """

    def __init__(
        self,
        store: tRNAseqDataStore,
        num_workers: int = 4,
        chunk_size: int = 1
    ):
        """
        Initialize batch processor.

        Args:
            store: tRNAseqDataStore instance
            num_workers: Number of parallel workers (default: 4)
            chunk_size: Batch size for processing (default: 1)
        """
        self.store = store
        self.num_workers = num_workers
        self.chunk_size = chunk_size
        self.logger = logging.getLogger("BatchProcessor")

    def process_batch(
        self,
        samples: List[str],
        func: Callable[[str], pd.DataFrame],
        output_type: str = 'charge',
        merge_results: bool = True,
        save_individual: bool = True,
        **kwargs
    ) -> Union[pd.DataFrame, List[pd.DataFrame]]:
        """
        Process multiple samples in parallel.

        Args:
            samples: List of sample IDs to process
            func: Processing function that takes sample_id and returns DataFrame
            output_type: Type of output ('charge', 'stats', or 'custom')
            merge_results: Whether to merge results into single DataFrame
            save_individual: Whether to save individual sample results
            **kwargs: Additional arguments passed to func

        Returns:
            Merged DataFrame if merge_results=True, else list of DataFrames

        Example:
            >>> def process_sample(sample_id, **kwargs):
            ...     # Load and process sample
            ...     return result_df
            >>>
            >>> results = processor.process_batch(
            ...     samples=['s1', 's2', 's3'],
            ...     func=process_sample
            ... )
        """
        tracker = ProgressTracker(total=len(samples), desc=f"Processing {output_type}")

        try:
            with Pool(processes=self.num_workers) as pool:
                results = []
                for sample_id in samples:
                    try:
                        result = func(sample_id, **kwargs)

                        if save_individual and output_type in ['charge', 'stats']:
                            if output_type == 'charge':
                                self.store.save_charge_data(result, sample_id=sample_id)
                            else:
                                self.store.save_stats(result, filename=f'stats_{sample_id}')

                        results.append(result)
                        tracker.update(f"Completed {sample_id}")

                    except Exception as e:
                        self.logger.error(f"Error processing {sample_id}: {str(e)}")
                        continue

            tracker.summary()

            if merge_results and results:
                return pd.concat(results, ignore_index=True)
            else:
                return results

        except Exception as e:
            self.logger.error(f"Batch processing failed: {str(e)}")
            raise

    def process_in_chunks(
        self,
        samples: List[str],
        func: Callable[[List[str]], pd.DataFrame],
        **kwargs
    ) -> pd.DataFrame:
        """
        Process samples in chunks (useful for memory efficiency).

        Args:
            samples: List of sample IDs
            func: Function that processes a list of samples
            **kwargs: Additional arguments passed to func

        Returns:
            Merged DataFrame of results

        Example:
            >>> def process_chunk(chunk_ids, store, **kwargs):
            ...     data = store.load_charge_data(samples=chunk_ids)
            ...     return analyze(data)
            >>>
            >>> results = processor.process_in_chunks(
            ...     samples=all_samples,
            ...     func=process_chunk
            ... )
        """
        tracker = ProgressTracker(
            total=len(samples),
            desc=f"Processing {len(samples)} samples in chunks"
        )

        results = []
        for i in range(0, len(samples), self.chunk_size):
            chunk = samples[i:i + self.chunk_size]
            try:
                result = func(chunk, store=self.store, **kwargs)
                results.append(result)
                tracker.update(f"Processed chunk {i//self.chunk_size + 1}")
            except Exception as e:
                self.logger.error(f"Error processing chunk {i}-{i+len(chunk)}: {str(e)}")
                continue

        tracker.summary()

        if results:
            return pd.concat(results, ignore_index=True)
        else:
            raise RuntimeError("No results produced from batch processing")

    def validate_results(
        self,
        results: Union[pd.DataFrame, List[pd.DataFrame]],
        required_columns: Optional[List[str]] = None,
        min_rows: int = 0
    ) -> bool:
        """
        Validate batch processing results.

        Args:
            results: DataFrame or list of DataFrames
            required_columns: List of columns that must be present
            min_rows: Minimum number of rows required

        Returns:
            True if valid, raises ValueError otherwise
        """
        if isinstance(results, list):
            dfs_to_check = results
        else:
            dfs_to_check = [results]

        for df in dfs_to_check:
            if len(df) < min_rows:
                raise ValueError(
                    f"DataFrame has {len(df)} rows, minimum required: {min_rows}"
                )

            if required_columns:
                missing = set(required_columns) - set(df.columns)
                if missing:
                    raise ValueError(f"Missing columns: {missing}")

        return True


def process_samples_parallel(
    samples: List[str],
    func: Callable[[str], pd.DataFrame],
    num_workers: int = 4,
    **kwargs
) -> List[pd.DataFrame]:
    """
    Process multiple samples in parallel using multiprocessing.

    Simple parallel processing without result storage.

    Args:
        samples: List of sample IDs
        func: Function that takes sample_id and returns DataFrame
        num_workers: Number of parallel workers
        **kwargs: Additional arguments passed to func

    Returns:
        List of DataFrames from each sample

    Example:
        >>> def analyze_sample(sample_id, **kwargs):
        ...     return pd.DataFrame({'data': [1, 2, 3]})
        >>>
        >>> results = process_samples_parallel(
        ...     samples=['s1', 's2', 's3'],
        ...     func=analyze_sample,
        ...     num_workers=4
        ... )
    """
    tracker = ProgressTracker(total=len(samples), desc="Parallel processing")
    results = []

    try:
        with Pool(processes=num_workers) as pool:
            for sample_id in samples:
                try:
                    result = func(sample_id, **kwargs)
                    results.append(result)
                    tracker.update(f"Completed {sample_id}")
                except Exception as e:
                    logging.error(f"Error processing {sample_id}: {str(e)}")
                    continue

        tracker.summary()
        return results

    except Exception as e:
        logging.error(f"Parallel processing failed: {str(e)}")
        raise


def merge_sample_results(
    results: List[pd.DataFrame],
    sample_ids: List[str],
    add_sample_column: bool = True,
    sample_col_name: str = 'sample_id'
) -> pd.DataFrame:
    """
    Merge results from multiple samples.

    Args:
        results: List of DataFrames
        sample_ids: List of sample IDs corresponding to results
        add_sample_column: Whether to add sample ID column
        sample_col_name: Name for sample ID column

    Returns:
        Merged DataFrame with optional sample ID column
    """
    if len(results) != len(sample_ids):
        raise ValueError(
            f"Number of results ({len(results)}) must match sample_ids ({len(sample_ids)})"
        )

    if add_sample_column:
        for result, sample_id in zip(results, sample_ids):
            result[sample_col_name] = sample_id

    merged = pd.concat(results, ignore_index=True)

    # Move sample column to front if added
    if add_sample_column and sample_col_name in merged.columns:
        cols = [sample_col_name] + [c for c in merged.columns if c != sample_col_name]
        merged = merged[cols]

    return merged
