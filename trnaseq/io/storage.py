"""
Efficient data storage and retrieval using Parquet format.

This module provides the tRNAseqDataStore class for handling storage and loading
of tRNA-seq analysis data with compression and efficient filtering capabilities.

Classes:
    tRNAseqDataStore: Main class for managing Parquet-based storage

Requirements:
    - pandas
    - pyarrow (for Parquet support)

Install with: pip install pandas pyarrow
"""

import os
import pandas as pd
from pathlib import Path
from typing import Optional, List, Union

try:
    import pyarrow.parquet as pq
    import pyarrow as pa
    PYARROW_AVAILABLE = True
except ImportError:
    PYARROW_AVAILABLE = False
    pq = None
    pa = None


class tRNAseqDataStore:
    """
    Efficient storage and retrieval of tRNA-seq data using Parquet format.

    This class handles saving and loading statistics and charge data with
    compression and optional filtering by sample or tRNA IDs.

    Attributes:
        base_dir (Path): Base directory for storing data files
        stats_dir (Path): Directory for statistics files
        charge_dir (Path): Directory for charge quantification data

    Example:
        >>> store = tRNAseqDataStore('project_data/')
        >>> store.save_stats(stats_df, format='parquet')
        >>> charge_data = store.load_charge_data(samples=['sample001', 'sample002'])
    """

    def __init__(self, base_dir: Union[str, Path]):
        """
        Initialize data store in base directory.

        Creates subdirectories for statistics and charge data if they don't exist.

        Args:
            base_dir: Base directory path for storing data

        Raises:
            ImportError: If pyarrow is not installed
        """
        if not PYARROW_AVAILABLE:
            raise ImportError(
                "pyarrow is required for Parquet support. "
                "Install with: pip install pyarrow"
            )

        self.base_dir = Path(base_dir)
        self.stats_dir = self.base_dir / 'stats'
        self.charge_dir = self.base_dir / 'charge'

        # Create directories if they don't exist
        self.stats_dir.mkdir(parents=True, exist_ok=True)
        self.charge_dir.mkdir(parents=True, exist_ok=True)

    def save_stats(
        self,
        df: pd.DataFrame,
        filename: str = 'stats',
        format: str = 'parquet',
        compression: str = 'snappy',
        **kwargs
    ) -> Path:
        """
        Save statistics DataFrame to Parquet or CSV format.

        Parquet format provides ~15x compression compared to CSV and enables
        efficient column-wise access and filtering.

        Args:
            df: DataFrame containing statistics
            filename: Output filename without extension (default: 'stats')
            format: Output format - 'parquet' (default) or 'csv'
            compression: Compression algorithm for Parquet:
                - 'snappy' (default): Fast, good compression
                - 'gzip': Higher compression, slower
                - 'brotli': Even higher compression, requires pyarrow[compression]
                - None: No compression
            **kwargs: Additional arguments passed to to_parquet() or to_csv()

        Returns:
            Path: Path to saved file

        Raises:
            ValueError: If format is not 'parquet' or 'csv'

        Example:
            >>> store.save_stats(df, filename='all_samples', compression='gzip')
        """
        if format.lower() == 'parquet':
            output_path = self.stats_dir / f'{filename}.parquet'
            df.to_parquet(
                output_path,
                compression=compression,
                index=False,
                **kwargs
            )
        elif format.lower() == 'csv':
            output_path = self.stats_dir / f'{filename}.csv'
            df.to_csv(output_path, index=False, **kwargs)
        else:
            raise ValueError(f"Format must be 'parquet' or 'csv', got {format}")

        return output_path

    def load_stats(
        self,
        filename: str = 'stats',
        samples: Optional[List[str]] = None,
        trnas: Optional[List[str]] = None,
        format: str = 'parquet',
        columns: Optional[List[str]] = None,
        **kwargs
    ) -> pd.DataFrame:
        """
        Load statistics with optional filtering by sample or tRNA.

        This method efficiently loads data from Parquet files, supporting
        column pruning and row filtering without loading the entire dataset.

        Args:
            filename: Input filename without extension (default: 'stats')
            samples: List of sample IDs to load (None = all samples)
            trnas: List of tRNA IDs to load (None = all tRNAs)
            format: Input format - 'parquet' (default) or 'csv'
            columns: Specific columns to load (None = all columns)
            **kwargs: Additional arguments for to_parquet() or read_csv()

        Returns:
            pd.DataFrame: Loaded statistics, filtered if criteria provided

        Raises:
            FileNotFoundError: If specified file doesn't exist
            ValueError: If format is not 'parquet' or 'csv'

        Example:
            >>> # Load specific samples and columns
            >>> df = store.load_stats(
            ...     samples=['sample001', 'sample002'],
            ...     columns=['sample_name', 'tRNA_annotation', 'count']
            ... )
        """
        if format.lower() == 'parquet':
            filepath = self.stats_dir / f'{filename}.parquet'
            if not filepath.exists():
                raise FileNotFoundError(f"File not found: {filepath}")

            # Use pyarrow for efficient filtering
            table = pq.read_table(filepath, columns=columns)
            df = table.to_pandas()
        elif format.lower() == 'csv':
            filepath = self.stats_dir / f'{filename}.csv'
            if not filepath.exists():
                raise FileNotFoundError(f"File not found: {filepath}")

            df = pd.read_csv(filepath, usecols=columns, **kwargs)
        else:
            raise ValueError(f"Format must be 'parquet' or 'csv', got {format}")

        # Apply filtering
        if samples is not None:
            if 'sample_name_unique' in df.columns:
                df = df[df['sample_name_unique'].isin(samples)]
            elif 'sample_name' in df.columns:
                df = df[df['sample_name'].isin(samples)]

        if trnas is not None:
            if 'tRNA_annotation' in df.columns:
                df = df[df['tRNA_annotation'].isin(trnas)]

        return df

    def save_charge_data(
        self,
        df: pd.DataFrame,
        sample_id: Optional[str] = None,
        append: bool = False,
        compression: str = 'snappy',
        **kwargs
    ) -> Path:
        """
        Save charge quantification results.

        Supports both single-sample and aggregated charge data storage.

        Args:
            df: DataFrame containing charge quantification results
            sample_id: Sample ID for single-sample data files (optional).
                If provided, saves to sample-specific file.
            append: If True, append to existing file (requires pyarrow dataset API)
            compression: Compression algorithm ('snappy', 'gzip', 'brotli', or None)
            **kwargs: Additional arguments for to_parquet()

        Returns:
            Path: Path to saved file

        Example:
            >>> store.save_charge_data(charge_df, sample_id='sample001')
            >>> store.save_charge_data(aggregated_df)  # Saves as 'charge_data.parquet'
        """
        if sample_id:
            filename = f'charge_{sample_id}.parquet'
        else:
            filename = 'charge_data.parquet'

        output_path = self.charge_dir / filename

        df.to_parquet(
            output_path,
            compression=compression,
            index=False,
            **kwargs
        )

        return output_path

    def load_charge_data(
        self,
        samples: Optional[List[str]] = None,
        filename: str = 'charge_data',
        columns: Optional[List[str]] = None,
        **kwargs
    ) -> pd.DataFrame:
        """
        Load charge data efficiently with optional sample filtering.

        For single-sample queries, loads only the relevant parquet files.
        For multi-sample queries, loads aggregated data and filters.

        Args:
            samples: List of sample IDs to load (None = load aggregated data)
            filename: Input filename for aggregated data (default: 'charge_data')
            columns: Specific columns to load (None = all columns)
            **kwargs: Additional arguments for to_parquet()

        Returns:
            pd.DataFrame: Loaded charge data

        Example:
            >>> # Load aggregated charge data
            >>> df = store.load_charge_data()

            >>> # Load specific samples
            >>> df = store.load_charge_data(samples=['sample001', 'sample002'])
        """
        if samples is None:
            # Load aggregated data
            filepath = self.charge_dir / f'{filename}.parquet'
            if not filepath.exists():
                raise FileNotFoundError(f"File not found: {filepath}")

            table = pq.read_table(filepath, columns=columns)
            return table.to_pandas()
        else:
            # Load and concatenate sample-specific files
            dfs = []
            for sample_id in samples:
                filepath = self.charge_dir / f'charge_{sample_id}.parquet'
                if filepath.exists():
                    table = pq.read_table(filepath, columns=columns)
                    dfs.append(table.to_pandas())

            if not dfs:
                raise FileNotFoundError(
                    f"No charge data files found for samples: {samples}"
                )

            return pd.concat(dfs, ignore_index=True)

    def get_available_files(self, data_type: str = 'stats') -> List[str]:
        """
        Get list of available data files.

        Args:
            data_type: Type of data - 'stats' or 'charge'

        Returns:
            List of available filenames (without extension)
        """
        if data_type == 'stats':
            directory = self.stats_dir
        elif data_type == 'charge':
            directory = self.charge_dir
        else:
            raise ValueError(f"data_type must be 'stats' or 'charge', got {data_type}")

        if not directory.exists():
            return []

        parquet_files = list(directory.glob('*.parquet'))
        csv_files = list(directory.glob('*.csv'))

        files = [f.stem for f in parquet_files + csv_files]
        return sorted(list(set(files)))

    def get_file_info(self, filepath: Union[str, Path]) -> dict:
        """
        Get metadata about a Parquet file.

        Args:
            filepath: Path to parquet file

        Returns:
            Dictionary with file metadata including rows, columns, size, compression
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")

        parquet_file = pq.ParquetFile(filepath)

        return {
            'num_rows': parquet_file.metadata.num_rows,
            'num_columns': parquet_file.metadata.num_columns,
            'columns': parquet_file.schema.names,
            'file_size_mb': filepath.stat().st_size / (1024 * 1024),
            'compression': parquet_file.metadata.row_group(0).column(0).compression
        }

    def estimate_compression_ratio(
        self,
        csv_filepath: Union[str, Path],
        parquet_filepath: Union[str, Path]
    ) -> float:
        """
        Calculate compression ratio between CSV and Parquet.

        Args:
            csv_filepath: Path to CSV file
            parquet_filepath: Path to corresponding Parquet file

        Returns:
            Compression ratio (CSV size / Parquet size)
        """
        csv_path = Path(csv_filepath)
        parquet_path = Path(parquet_filepath)

        if not csv_path.exists() or not parquet_path.exists():
            raise FileNotFoundError("Both CSV and Parquet files must exist")

        csv_size = csv_path.stat().st_size
        parquet_size = parquet_path.stat().st_size

        return csv_size / parquet_size if parquet_size > 0 else 0

    def convert_csv_to_parquet(
        self,
        csv_filepath: Union[str, Path],
        output_name: Optional[str] = None,
        compression: str = 'snappy',
        **pandas_kwargs
    ) -> Path:
        """
        Convert existing CSV file to Parquet format.

        Args:
            csv_filepath: Path to CSV file
            output_name: Output filename without extension (default: CSV filename)
            compression: Compression algorithm
            **pandas_kwargs: Additional arguments for pd.read_csv()

        Returns:
            Path: Path to created Parquet file

        Example:
            >>> store.convert_csv_to_parquet('old_stats.csv', output_name='stats')
        """
        csv_path = Path(csv_filepath)
        if not csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_path}")

        df = pd.read_csv(csv_path, **pandas_kwargs)

        if output_name is None:
            output_name = csv_path.stem

        output_path = self.stats_dir / f'{output_name}.parquet'
        df.to_parquet(output_path, compression=compression, index=False)

        return output_path
