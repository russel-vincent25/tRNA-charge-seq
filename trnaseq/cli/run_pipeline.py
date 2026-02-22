#!/usr/bin/env python3
"""
Unified Preprocessing Pipeline (Stages 0-5)
============================================

Combines all preprocessing and analysis steps into a single script:
- Stage 0a: Adapter Removal & Merging
- Stage 0b: Barcode Splitting
- Stage 0c: UMI Trimming
- Stage 1: SWIPE Alignment
- Stage 2: Stats Collection
- Stage 3: Charge Quantification (optional)
- Stage 4: Parquet Storage (optional)
- Stage 5: QC Summary Report (optional)

Output:
- inp_file_df.xlsx (input file summary)
- sample_df.xlsx (sample information with QC metrics)
- ALL_stats_aggregate.csv (combined statistics for charge quantification)
- charge_analysis/ (charge quantification results, if enabled)
- parquet_data/ (parquet-format data, if enabled)
- QC_summary.csv, QC_report.html (QC dashboard, if enabled)

Usage:
    python -m trnaseq pipeline \
        --config config.yaml \
        --project-dir /path/to/project/ \
        --n-jobs 8 \
        --parquet

Author: Based on projects/example/process_data.ipynb
Date: 2026-02-11
"""

import os
import sys
import time
import json
from pathlib import Path
from datetime import datetime
import argparse
import yaml
import pandas as pd
import numpy as np

# Add repo root to path for src imports
repo_path = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_path))

from src.misc import (
    index_to_sample_df,
    downsample_raw_input,
    read_tRNAdb_info,
    sample_df_to_dict
)
from src.read_processing import AR_merge, BC_split, UMI_trim
from src.alignment import SWIPE_align
from src.stats_collection import STATS_collection

# Optional imports for extended analysis stages
try:
    from trnaseq.charge.quantifier import ChargeQuantifier
    CHARGE_AVAILABLE = True
except ImportError:
    CHARGE_AVAILABLE = False

try:
    from trnaseq.io.storage import tRNAseqDataStore, PYARROW_AVAILABLE
    STORAGE_AVAILABLE = PYARROW_AVAILABLE
except ImportError:
    STORAGE_AVAILABLE = False

try:
    from trnaseq.qc.report import QCReportGenerator
    QC_AVAILABLE = True
except ImportError:
    QC_AVAILABLE = False

try:
    from trnaseq.modifications.positional import PositionalExtractor
    from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer
    from trnaseq.modifications.modification_caller import ModificationCaller
    from trnaseq.modifications.modomics import MODOMICSAnnotator
    MODIFICATIONS_AVAILABLE = True
except ImportError:
    MODIFICATIONS_AVAILABLE = False


# ------------------------------------------------------------------
# Path fields in config that should be resolved relative to project_dir
# ------------------------------------------------------------------
_PATH_FIELDS = ['sample_list', 'index_list', 'SWIPE_score_mat', 'common_seqs', 'adapter_sequences']
_DICT_PATH_FIELDS = ['tRNA_database']  # dict values are paths


class PreprocessingPipeline:
    """
    Unified preprocessing pipeline for tRNA-charge-seq

    Runs stages 0a-5 and outputs key files for downstream analysis.
    """

    def __init__(self, config_file, project_dir, n_jobs=4, sample_index=None,
                 threads_per_job=None):
        """
        Initialize pipeline

        Parameters:
            config_file: YAML configuration file
            project_dir: Project directory (must contain data/{seq_dir}/ for stages 0-1)
            n_jobs: Number of parallel jobs
            sample_index: Process only sample at this 0-based index (for SLURM array jobs)
            threads_per_job: Threads per subprocess (AR/SWIPE). Default: from config or 2.
        """
        self.config_file = Path(config_file).resolve()
        self.project_dir = Path(project_dir).resolve()
        self.n_jobs = n_jobs
        self.sample_index = sample_index
        self.log_file = self.project_dir / "preprocessing.log"

        # Create output directory
        self.project_dir.mkdir(parents=True, exist_ok=True)

        # Load configuration
        with open(self.config_file, 'r') as f:
            self.config = yaml.safe_load(f)

        # Resolve relative paths in config against project_dir
        self._resolve_config_paths()

        # Resolve threads_per_job: CLI arg > config > default 2
        if threads_per_job is not None:
            self.threads_per_job = threads_per_job
        else:
            self.threads_per_job = self.config.get('threads_per_job', 2)

        # Initialize status tracking
        self.status = {
            'start_time': datetime.now(),
            'stages_completed': []
        }

        # Computing metrics
        self.stage_timings = {}
        self.file_metrics = {}

    def _resolve_config_paths(self):
        """Resolve relative paths in config against project_dir.

        A path is treated as relative if it doesn't start with '/'.
        Absolute paths pass through unchanged. None/null values are skipped.
        """
        for key in _PATH_FIELDS:
            val = self.config.get(key)
            if val is not None and not os.path.isabs(val):
                self.config[key] = str(self.project_dir / val)

        for key in _DICT_PATH_FIELDS:
            d = self.config.get(key)
            if isinstance(d, dict):
                for k, v in d.items():
                    if v is not None and not os.path.isabs(v):
                        d[k] = str(self.project_dir / v)

    def log(self, message, level="INFO"):
        """Log message to file and console"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_msg = f"[{timestamp}] {level}: {message}"

        # Print to console
        print(log_msg)

        # Write to log file
        with open(self.log_file, 'a') as f:
            f.write(log_msg + '\n')

    def load_sample_info(self):
        """Load sample list and index information"""
        self.log("Loading sample information...")

        # Load sample list
        sample_list = self.config['sample_list']
        self.sample_df = pd.read_excel(sample_list)

        # Load index list
        index_list = self.config['index_list']
        index_df = pd.read_excel(index_list)

        # Add barcode and index sequences to sample_df
        self.sample_df = index_to_sample_df(self.sample_df, index_df)

        # Create input file DataFrame
        self.inp_file_df = self.sample_df[[
            'fastq_mate1_filename', 'fastq_mate2_filename',
            'P5_index', 'P7_index', 'P5_index_seq', 'P7_index_seq'
        ]].drop_duplicates().reset_index(drop=True)

        self.log(f"Loaded {len(self.sample_df)} samples")
        self.log(f"Input files: {len(self.inp_file_df)}")

    def setup_directories(self):
        """Setup directory structure"""
        self.log("Setting up directories...")

        self.dir_dict = {
            'NBdir': str(self.project_dir),
            'data_dir': 'data',
            'seq_dir': self.config.get('seq_dir', 'raw_fastq'),
            'AdapterRemoval_dir': 'AdapterRemoval',
            'BC_dir': 'BC_split',
            'UMI_dir': 'UMI_trimmed',
            'align_dir': 'SWalign',
            'stats_dir': 'stats_collection',
            'charge_dir': 'charge_analysis',
            'parquet_dir': 'parquet_data',
        }

        # Create data directory
        (self.project_dir / 'data').mkdir(exist_ok=True)

        # Create subdirectories under data/
        for dir_name in ['AdapterRemoval', 'BC_split', 'UMI_trimmed',
                         'SWalign', 'stats_collection']:
            (self.project_dir / 'data' / dir_name).mkdir(exist_ok=True)

        # Create analysis directories
        for dir_name in ['charge_analysis', 'parquet_data']:
            (self.project_dir / dir_name).mkdir(exist_ok=True)

    def _get_working_sample_df(self):
        """Return sample_df subset if --sample-index is set, else full."""
        if self.sample_index is not None:
            return self.sample_df.iloc[[self.sample_index]].copy()
        return self.sample_df

    def stage_0a_adapter_removal(self):
        """Stage 0a: Adapter Removal & Merging"""
        self.log("=" * 60)
        self.log("Stage 0a: Adapter Removal & Merging")
        self.log("=" * 60)

        MIN_READ_LEN = self.config.get('min_read_len', 39)

        # Run AdapterRemoval
        AR_obj = AR_merge(
            self.dir_dict, self.inp_file_df, MIN_READ_LEN,
            AR_threads=self.threads_per_job,
            overwrite_dir=self.config.get('overwrite', True)
        )
        self.inp_file_df = AR_obj.run_parallel(
            n_jobs=self.n_jobs,
            overwrite=self.config.get('overwrite', True)
        )

        # Log statistics
        for _, row in self.inp_file_df.iterrows():
            self.log(f"  {row['fastq_mate1_filename']}: "
                    f"{row['percent_successfully_merged']:.1f}% merged")

        # --- Merge QC probe: check UMI + barcode detectability ---
        try:
            barcodes = list(self.sample_df['barcode_seq'].unique())
            umi_mode = self.config.get('umi_trim_mode', 'anchored')
            probe_kwargs = {'umi_mode': umi_mode}

            if umi_mode == 'anchored':
                adapter_seqs_path = self.config.get('adapter_sequences', None)
                if adapter_seqs_path is None:
                    adapter_seqs_path = str(self.project_dir / 'utils' / 'adapter_sequences.yaml')
                with open(adapter_seqs_path) as f:
                    adapter_seqs = yaml.safe_load(f)
                anchor_name = self.config.get('umi_anchor')
                if anchor_name and anchor_name in adapter_seqs:
                    probe_kwargs['anchor_seq'] = adapter_seqs[anchor_name]
                    probe_kwargs['max_stagger'] = self.config.get('umi_max_stagger', 3)
                    probe_kwargs['anchor_max_dist'] = self.config.get('umi_anchor_max_dist', 1)

            probe_results = AR_obj.qc_probe(barcodes, **probe_kwargs)

            self.log(f"  Merge QC probe ({self.config.get('qc_probe_reads', 10000)} reads per file pair):")
            total_umi = total_bc = total_both = total_n = 0
            for basename, stats in probe_results.items():
                self.log(f"    {basename}: UMI={stats['pct_umi']}% "
                         f"BC={stats['pct_bc']}% both={stats['pct_both']}%")
                if 'stagger_counts' in stats:
                    self.log(f"    {basename}: stagger distribution: {stats['stagger_counts']}")
                total_umi += stats['pct_umi'] * stats['n_probed']
                total_bc += stats['pct_bc'] * stats['n_probed']
                total_both += stats['pct_both'] * stats['n_probed']
                total_n += stats['n_probed']

            if total_n > 0:
                avg_umi = total_umi / total_n
                avg_bc = total_bc / total_n
                avg_both = total_both / total_n
                self.log(f"  Average: UMI={avg_umi:.1f}% BC={avg_bc:.1f}% both={avg_both:.1f}%")
                if avg_umi < 50:
                    self.log("  WARNING: Low UMI detection rate (<50%). "
                             "Check umi_trim_mode / umi_anchor config.", level="WARN")
                if avg_bc < 50:
                    self.log("  WARNING: Low barcode detection rate (<50%). "
                             "Check barcode sequences in sample list.", level="WARN")
        except Exception as e:
            self.log(f"  Merge QC probe skipped: {e}", level="WARN")

        self.status['stages_completed'].append('0a')

    def stage_0b_barcode_split(self):
        """Stage 0b: Barcode Splitting"""
        self.log("=" * 60)
        self.log("Stage 0b: Barcode Splitting")
        self.log("=" * 60)

        # Run barcode splitting
        BCsplit_obj = BC_split(
            self.dir_dict, self.sample_df, self.inp_file_df,
            overwrite_dir=self.config.get('overwrite', True)
        )
        self.sample_df, self.inp_file_df = BCsplit_obj.run_parallel(
            n_jobs=self.n_jobs
        )

        # Log statistics
        for _, row in self.inp_file_df.iterrows():
            if row['percent_BC-mapped'] < 80:
                self.log(f"  WARNING: Low barcode mapping: "
                        f"{row['percent_BC-mapped']:.1f}%", level="WARN")
            else:
                self.log(f"  Barcode mapping: {row['percent_BC-mapped']:.1f}%")

        self.status['stages_completed'].append('0b')

    def stage_0c_umi_trim(self):
        """Stage 0c: UMI Trimming"""
        self.log("=" * 60)
        self.log("Stage 0c: UMI Trimming")
        self.log("=" * 60)

        # Run UMI trimming
        downsample_percentile = self.config.get('downsample_percentile', None)
        downsample_absolute = self.config.get('downsample_absolute', None)

        # UMI trim mode: 'anchored' (default) or 'pyrimidine' (legacy)
        umi_mode = self.config.get('umi_trim_mode', 'anchored')
        anchor_seq = None
        if umi_mode == 'anchored':
            adapter_seqs_path = self.config.get('adapter_sequences', None)
            if adapter_seqs_path is None:
                adapter_seqs_path = str(self.project_dir / 'utils' / 'adapter_sequences.yaml')
            with open(adapter_seqs_path) as f:
                adapter_seqs = yaml.safe_load(f)
            anchor_name = self.config.get('umi_anchor')
            if anchor_name is None:
                raise ValueError("Config field 'umi_anchor' is required when umi_trim_mode='anchored'")
            if anchor_name not in adapter_seqs:
                raise ValueError(f"Anchor '{anchor_name}' not found in {adapter_seqs_path}. "
                                 f"Available: {list(adapter_seqs.keys())}")
            anchor_seq = adapter_seqs[anchor_name]
            self.log(f"  UMI mode: anchored (anchor={anchor_name}: {anchor_seq})")
        else:
            self.log(f"  UMI mode: pyrimidine (legacy)")

        working_df = self._get_working_sample_df()
        UMItrim_obj = UMI_trim(
            self.dir_dict, working_df,
            mode=umi_mode,
            anchor_seq=anchor_seq,
            max_stagger=self.config.get('umi_max_stagger', 3),
            anchor_max_dist=self.config.get('umi_anchor_max_dist', 1),
            overwrite_dir=self.config.get('overwrite', True) if self.sample_index is None else False,
            downsample_percentile=downsample_percentile,
            downsample_absolute=downsample_absolute
        )
        result_df = UMItrim_obj.run_parallel(n_jobs=self.n_jobs)
        if self.sample_index is not None:
            # Merge single-sample results back into full sample_df
            merge_cols = [c for c in result_df.columns if c not in self.sample_df.columns or c == 'sample_name_unique']
            extra_cols = [c for c in result_df.columns if c not in self.sample_df.columns]
            if extra_cols:
                self.sample_df = self.sample_df.merge(
                    result_df[['sample_name_unique'] + extra_cols],
                    on='sample_name_unique', how='left')
        else:
            self.sample_df = result_df

        # Log statistics
        avg_umi_pct = self.sample_df['percent_UMI_obs-vs-exp'].mean()
        self.log(f"  Average UMI obs/exp: {avg_umi_pct:.1f}%")

        if avg_umi_pct < 80:
            self.log(f"  WARNING: Low UMI diversity suggests library prep bottleneck",
                    level="WARN")

        self.status['stages_completed'].append('0c')

    def stage_1_alignment(self):
        """Stage 1: SWIPE Alignment"""
        self.log("=" * 60)
        self.log("Stage 1: SWIPE Alignment")
        self.log("=" * 60)

        # Load tRNA database info
        tRNA_database = self.config['tRNA_database']
        self.tRNA_data = read_tRNAdb_info(tRNA_database)

        # Setup alignment parameters
        SWIPE_score_mat = self.config['SWIPE_score_mat']
        common_seqs = self.config.get('common_seqs', None)
        MIN_SCORE_ALIGN = self.config.get('min_score_align', 15)

        # Run alignment
        working_df = self._get_working_sample_df()
        align_obj = SWIPE_align(
            self.dir_dict, tRNA_database, working_df,
            SWIPE_score_mat,
            gap_penalty=self.config.get('gap_penalty', 6),
            extension_penalty=self.config.get('extension_penalty', 3),
            min_score_align=MIN_SCORE_ALIGN,
            common_seqs=common_seqs,
            overwrite_dir=self.config.get('overwrite', True) if self.sample_index is None else False,
            SWIPE_threads=self.threads_per_job,
        )
        result_df = align_obj.run_parallel(n_jobs=self.n_jobs)
        if self.sample_index is not None:
            merge_cols = [c for c in result_df.columns if c not in self.sample_df.columns]
            if merge_cols:
                self.sample_df = self.sample_df.merge(
                    result_df[['sample_name_unique'] + merge_cols],
                    on='sample_name_unique', how='left')
        else:
            self.sample_df = result_df

        # Log statistics
        avg_mapping = self.sample_df['Mapping_percent'].mean()
        self.log(f"  Average mapping rate: {avg_mapping:.1f}%")

        if avg_mapping < 80:
            self.log(f"  WARNING: Low mapping rate", level="WARN")

        self.status['stages_completed'].append('1')

    def stage_2_stats_collection(self):
        """Stage 2: Stats Collection"""
        self.log("=" * 60)
        self.log("Stage 2: Stats Collection")
        self.log("=" * 60)

        common_seqs = self.config.get('common_seqs', None)

        # Check if common-seq-obs files actually exist in the SWalign directory
        # (common_seqs config points to the fasta used during alignment,
        #  but STATS_collection needs the *result* files to exist)
        if common_seqs is not None:
            align_dir = self.project_dir / 'data' / self.dir_dict['align_dir']
            common_obs_files = list(align_dir.glob('*_common-seq-obs.json*'))
            if not common_obs_files:
                self.log("  WARNING: common_seqs configured but no common-seq-obs files "
                        "found in SWalign/. Setting common_seqs=None.", level="WARN")
                common_seqs = None

        # Run stats collection
        stats_obj = STATS_collection(
            self.dir_dict, self.tRNA_data, self.sample_df,
            common_seqs=common_seqs,
            overwrite_dir=self.config.get('overwrite', True)
        )
        self.stats_df = stats_obj.run_parallel(n_jobs=self.n_jobs)

        # Copy fragment counts to sample_df
        self.fragment_counts = getattr(stats_obj, 'fragment_counts', {})
        if self.fragment_counts:
            for col in ['N_full_length', 'N_rt_dropoff', 'N_5p_fragment', 'N_degraded', 'N_total_aligned']:
                self.sample_df[col] = self.sample_df['sample_name_unique'].map(
                    lambda sn, c=col: self.fragment_counts.get(sn, {}).get(c, 0)
                )

            # Recalculate Mapping_percent from fragment counts:
            # N_total_aligned / N_after_trim * 100
            if 'N_after_trim' in self.sample_df.columns:
                self.sample_df['Mapping_percent'] = np.where(
                    self.sample_df['N_after_trim'] > 0,
                    self.sample_df['N_total_aligned'] / self.sample_df['N_after_trim'] * 100,
                    0.0
                )

        self.log(f"  Collected stats for {len(self.stats_df)} entries")

        self.status['stages_completed'].append('2')

    def stage_3_charge_quantification(self):
        """Stage 3: Charge Quantification"""
        self.log("=" * 60)
        self.log("Stage 3: Charge Quantification")
        self.log("=" * 60)

        if not CHARGE_AVAILABLE:
            self.log("WARNING: ChargeQuantifier not available. Skipping stage 3.", level="WARN")
            self.log("  Install with: pip install -e . (from repo root)", level="WARN")
            return

        # Get stats file
        stats_file = self.project_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
        if not stats_file.exists():
            self.log(f"ERROR: Stats file not found: {stats_file}", level="ERROR")
            return

        # Get charge quantification settings
        charge_count = self.config.get('charge_count', 'count')
        charge_levels = self.config.get('charge_levels', ['transcript', 'codon', 'aa'])
        include_synthetic = self.config.get('include_synthetic', False)
        include_mito = self.config.get('include_mito', True)

        self.log(f"  Input: {stats_file.name}")
        self.log(f"  Charge count column: {charge_count}")
        self.log(f"  Levels: {', '.join(charge_levels)}")

        try:
            # Initialize quantifier
            quantifier = ChargeQuantifier(
                stats_csv=str(stats_file),
                charge_count=charge_count,
                RPM_count=charge_count
            )

            self.log(f"  Loaded {len(quantifier.stats_df)} alignment records")

            # Create output directory
            charge_dir = self.project_dir / 'charge_analysis'
            charge_dir.mkdir(exist_ok=True)

            # Quantify and export for each level
            self.charge_results = {}
            for level in charge_levels:
                self.log(f"  Quantifying at {level} level...")

                df = quantifier.quantify_all(
                    level=level,
                    include_synthetic=include_synthetic,
                    include_mito=include_mito
                )

                # Save to CSV
                output_file = charge_dir / f'charge_df_{level}.csv'
                df.to_csv(output_file, index=False)
                self.charge_results[level] = df

                self.log(f"    Saved: {output_file.name} ({len(df)} entries)")

            # Generate summary statistics
            self.log("  Generating summary statistics...")
            summary_list = []
            for level in charge_levels:
                summary = quantifier.get_summary_statistics(level=level)
                summary['level'] = level
                summary_list.append(summary)

            if summary_list:
                summary_df = pd.concat(summary_list, ignore_index=True)
                summary_file = charge_dir / 'charge_summary.csv'
                summary_df.to_csv(summary_file, index=False)
                self.charge_summary_df = summary_df
                self.log(f"    Saved: {summary_file.name}")

            self.status['stages_completed'].append('3')
            self.log("  Charge quantification complete!")

        except Exception as e:
            self.log(f"ERROR in charge quantification: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def stage_4_parquet_storage(self):
        """Stage 4: Parquet Storage (optional)"""
        self.log("=" * 60)
        self.log("Stage 4: Parquet Storage")
        self.log("=" * 60)

        if not STORAGE_AVAILABLE:
            self.log("WARNING: tRNAseqDataStore not available. Skipping stage 4.", level="WARN")
            self.log("  Install with: pip install pyarrow", level="WARN")
            return

        try:
            # Initialize data store
            parquet_dir = self.project_dir / 'parquet_data'
            store = tRNAseqDataStore(str(parquet_dir))

            compression = self.config.get('parquet_compression', 'snappy')
            self.log(f"  Compression: {compression}")

            # Convert stats CSV to Parquet
            stats_file = self.project_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
            if stats_file.exists():
                self.log(f"  Converting stats to Parquet...")
                stats_df = pd.read_csv(stats_file)
                parquet_file = store.save_stats(
                    stats_df,
                    filename='ALL_stats_aggregate',
                    format='parquet',
                    compression=compression
                )

                # Calculate compression ratio
                csv_size = stats_file.stat().st_size / (1024 * 1024)
                parquet_size = parquet_file.stat().st_size / (1024 * 1024)
                ratio = csv_size / parquet_size if parquet_size > 0 else 0

                self.log(f"    CSV: {csv_size:.1f} MB")
                self.log(f"    Parquet: {parquet_size:.1f} MB")
                self.log(f"    Compression ratio: {ratio:.1f}x")

            # Convert charge CSVs to Parquet if they exist
            charge_dir = self.project_dir / 'charge_analysis'
            if charge_dir.exists():
                charge_parquet_dir = parquet_dir / 'charge'
                charge_parquet_dir.mkdir(parents=True, exist_ok=True)
                for csv_file in charge_dir.glob('charge_df_*.csv'):
                    self.log(f"  Converting {csv_file.name} to Parquet...")
                    df = pd.read_csv(csv_file)
                    out_path = charge_parquet_dir / csv_file.with_suffix('.parquet').name
                    df.to_parquet(out_path, compression=compression)
                    self.log(f"    Saved: {out_path.name}")

            self.status['stages_completed'].append('4')
            self.log("  Parquet conversion complete!")

        except Exception as e:
            self.log(f"ERROR in Parquet storage: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def stage_5_qc_report(self):
        """Stage 5: QC Summary Report

        Generates QC_summary.csv and QC_report.html from data already
        computed by earlier stages. No new computation — just formatting.
        """
        self.log("=" * 60)
        self.log("Stage 5: QC Summary Report")
        self.log("=" * 60)

        if not QC_AVAILABLE:
            self.log("WARNING: QCReportGenerator not available. Skipping stage 5.", level="WARN")
            return

        try:
            # Load stats_df if not already in memory
            stats_df = getattr(self, 'stats_df', None)
            if stats_df is None:
                stats_file = (self.project_dir / 'data' / 'stats_collection'
                              / 'ALL_stats_aggregate.csv')
                if stats_file.exists():
                    stats_df = pd.read_csv(stats_file)

            # Load charge summary if not already in memory
            charge_summary = getattr(self, 'charge_summary_df', None)
            if charge_summary is None:
                cs_file = self.project_dir / 'charge_analysis' / 'charge_summary.csv'
                if cs_file.exists():
                    charge_summary = pd.read_csv(cs_file)

            # Ensure we have sample_df and inp_file_df
            sample_df = getattr(self, 'sample_df', None)
            inp_file_df = getattr(self, 'inp_file_df', pd.DataFrame())

            if sample_df is None:
                self.log("ERROR: sample_df not available for QC report", level="ERROR")
                return

            qc = QCReportGenerator(
                project_dir=self.project_dir,
                sample_df=sample_df,
                inp_file_df=inp_file_df,
                charge_summary_df=charge_summary,
                stats_df=stats_df,
                bc_dir=self.project_dir / 'data' / 'BC_split',
            )

            # QC_summary.csv
            summary_path = self.project_dir / 'QC_summary.csv'
            qc_summary = qc.save_summary_csv(summary_path)
            self.log(f"  Saved: {summary_path.name} ({len(qc_summary)} samples)")

            # QC_report.html
            report_path = self.project_dir / 'QC_report.html'
            qc.generate_html_report(report_path)
            self.log(f"  Saved: {report_path.name}")

            self.status['stages_completed'].append('5')
            self.log("  QC report generation complete!")

        except Exception as e:
            self.log(f"ERROR in QC report: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def stage_6_modification_analysis(self):
        """Stage 6: Modification Analysis (optional)

        Per-position RT signature extraction, MODOMICS annotation, and
        optional novel modification discovery.  Requires SWalign JSONs from
        stage 1.
        """
        self.log("=" * 60)
        self.log("Stage 6: Modification Analysis")
        self.log("=" * 60)

        if not MODIFICATIONS_AVAILABLE:
            self.log("WARNING: Modification analysis modules not available. "
                     "Skipping stage 6.", level="WARN")
            return

        try:
            # Determine reference FASTA (first species in tRNA_database config)
            tRNA_database = self.config['tRNA_database']
            ref_fasta = list(tRNA_database.values())[0]

            organism = self.config.get('organism', 'Escherichia coli')
            min_coverage = self.config.get('modification_min_coverage', 50)
            discover_novel = self.config.get('discover_novel_modifications', False)
            use_api = not self.config.get('no_modomics', True)

            json_dir = self.project_dir / 'data' / self.dir_dict['align_dir']
            output_dir = self.project_dir / 'modification_analysis'
            output_dir.mkdir(exist_ok=True)

            sample_names = self.sample_df['sample_name_unique'].tolist()

            self.log(f"  Reference: {ref_fasta}")
            self.log(f"  Organism: {organism}")
            self.log(f"  Samples: {len(sample_names)}")
            self.log(f"  Novel discovery: {discover_novel}")

            # Extract PSCMs
            extractor = PositionalExtractor(ref_fasta)
            all_pscm = extractor.run_parallel(
                json_dir, sample_names, n_jobs=self.n_jobs
            )

            # MODOMICS
            annotator = MODOMICSAnnotator(organism)
            mods_df = annotator.get_modifications(use_api=use_api)
            self.log(f"  Loaded {len(mods_df)} known modification entries")

            # Analyze each sample
            analyzer = RTSignatureAnalyzer(
                min_coverage=min_coverage, verbose=False
            )
            analyzer.load_reference(ref_fasta)
            caller = ModificationCaller(organism=organism)

            for sample_name, pscm_dict in all_pscm.items():
                pscm_dfs = analyzer.load_pscm_from_positional(pscm_dict)

                rt_profile = extractor.compute_rt_profile(pscm_dict)
                mm_profile = extractor.compute_mismatch_profile(pscm_dict)

                sample_dir = output_dir / sample_name
                sample_dir.mkdir(exist_ok=True)

                rt_profile.to_parquet(sample_dir / 'rt_profile.parquet', index=False)
                mm_profile.to_parquet(sample_dir / 'mismatch_profile.parquet', index=False)

                self.log(f"  {sample_name}: {len(pscm_dict)} tRNAs processed")

            self.status['stages_completed'].append('6')
            self.log("  Modification analysis complete!")

        except Exception as e:
            self.log(f"ERROR in modification analysis: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def _file_size_mb(self, path):
        """Return file size in MB, or None if file doesn't exist."""
        p = Path(path)
        if p.exists():
            return round(p.stat().st_size / (1024 * 1024), 3)
        return None

    def _collect_file_metrics(self):
        """Collect per-sample file sizes from pipeline output directories."""
        data_dir = self.project_dir / 'data'
        samples = self.sample_df['sample_name_unique'].tolist()

        for sample in samples:
            self.file_metrics[sample] = {
                'bc_split_fastq_MB': self._file_size_mb(
                    data_dir / 'BC_split' / f'{sample}.fastq.bz2'),
                'umi_trimmed_fastq_MB': self._file_size_mb(
                    data_dir / 'UMI_trimmed' / f'{sample}_UMI-trimmed.fastq.bz2'),
                'alignment_json_MB': self._file_size_mb(
                    data_dir / 'SWalign' / f'{sample}_SWalign.json.bz2'),
                'per_read_stats_MB': self._file_size_mb(
                    data_dir / 'stats_collection' / f'{sample}_stats.csv.bz2'),
            }

    def _save_metrics(self):
        """Save computing_metrics.csv (or per-sample JSON in single-sample mode)."""
        if self.sample_index is not None:
            # Single-sample SLURM mode: save per-sample JSON
            metrics_dir = self.project_dir / 'data' / 'metrics'
            metrics_dir.mkdir(exist_ok=True)

            sample_name = self.sample_df.iloc[self.sample_index]['sample_name_unique']
            sample_metrics = {
                'sample_name_unique': sample_name,
                **self.file_metrics.get(sample_name, {}),
                **{f'{k}_sec': round(v, 2) for k, v in self.stage_timings.items()},
            }

            out_path = metrics_dir / f'{sample_name}_metrics.json'
            with open(out_path, 'w') as f:
                json.dump(sample_metrics, f, indent=2)
            self.log(f"  Saved per-sample metrics: {out_path}")
        else:
            # Full pipeline or aggregation mode: save CSV
            self._aggregate_metrics_csv()

    def _aggregate_metrics_csv(self):
        """Aggregate per-sample JSONs (if any) and pipeline metrics into CSV."""
        metrics_dir = self.project_dir / 'data' / 'metrics'
        rows = []

        # Try loading per-sample JSONs from SLURM array jobs
        if metrics_dir.exists():
            for json_file in sorted(metrics_dir.glob('*_metrics.json')):
                with open(json_file) as f:
                    rows.append(json.load(f))

        if rows:
            # Merge with any file metrics collected in this run
            df = pd.DataFrame(rows)
        else:
            # No per-sample JSONs — build from current run's data
            self._collect_file_metrics()
            samples = self.sample_df['sample_name_unique'].tolist()
            for sample in samples:
                row = {'sample_name_unique': sample}
                row.update(self.file_metrics.get(sample, {}))
                rows.append(row)
            df = pd.DataFrame(rows)

        # Add pipeline-level timings as columns (same value for all rows)
        for stage, duration in self.stage_timings.items():
            col = f'{stage}_sec'
            if col not in df.columns:
                df[col] = round(duration, 2)

        # Add aggregate stats file size
        agg_stats = self.project_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
        if agg_stats.exists():
            df['aggregate_stats_MB'] = self._file_size_mb(agg_stats)

        out_path = self.project_dir / 'computing_metrics.csv'
        df.to_csv(out_path, index=False)
        self.log(f"  Saved computing metrics: {out_path} ({len(df)} samples)")

    def save_outputs(self):
        """Save output files"""
        self.log("=" * 60)
        self.log("Saving outputs")
        self.log("=" * 60)

        # Save inp_file_df
        inp_file_output = self.project_dir / 'inp_file_df.xlsx'
        self.inp_file_df.to_excel(inp_file_output, index=False)
        self.log(f"  Saved: {inp_file_output}")

        # Save sample_df
        sample_df_output = self.project_dir / 'sample_df.xlsx'
        self.sample_df.to_excel(sample_df_output, index=False)
        self.log(f"  Saved: {sample_df_output}")

        # ALL_stats_aggregate.csv should already be in stats_collection/
        stats_file = self.project_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
        if stats_file.exists():
            # Copy to output root for easy access
            import shutil
            shutil.copy(stats_file, self.project_dir / 'ALL_stats_aggregate.csv')
            self.log(f"  Saved: ALL_stats_aggregate.csv")
        else:
            self.log(f"  WARNING: ALL_stats_aggregate.csv not found", level="WARN")

        # Copy charge analysis results to output root
        if hasattr(self, 'charge_results') and self.charge_results:
            charge_dir = self.project_dir / 'charge_analysis'
            for level, df in self.charge_results.items():
                output_file = self.project_dir / f'charge_df_{level}.csv'
                import shutil
                shutil.copy(
                    charge_dir / f'charge_df_{level}.csv',
                    output_file
                )
                self.log(f"  Saved: {output_file.name}")

    @staticmethod
    def parse_stages(stages_str):
        """
        Parse --stages argument into a set of stage identifiers.

        Accepts:
            'all'       -> all stages
            '0-2'       -> stages 0a, 0b, 0c, 1, 2
            '3-5'       -> stages 3, 4, 5
            '0a,1,3'    -> specific stages
            '2,3,5'     -> specific stages

        Returns:
            Set of stage identifiers: {'0a', '0b', '0c', '1', '2', '3', '4', '5'}
        """
        ALL_STAGES = {'0a', '0b', '0c', '1', '2', '3', '4', '5', '6'}

        if stages_str == 'all':
            return ALL_STAGES

        stages = set()
        for part in stages_str.split(','):
            part = part.strip()
            if '-' in part and not part.startswith('0'):
                # Range like '3-5'
                start, end = part.split('-')
                # Expand numeric ranges
                for i in range(int(start), int(end) + 1):
                    stages.add(str(i))
            elif part == '0-2' or part == '0-1':
                # Special: includes sub-stages
                stages.update({'0a', '0b', '0c'})
                if part == '0-2':
                    stages.update({'1', '2'})
                else:
                    stages.add('1')
            elif part in ALL_STAGES:
                stages.add(part)
            elif part == '0':
                stages.update({'0a', '0b', '0c'})
            else:
                raise ValueError(f"Unknown stage: '{part}'. Valid: {sorted(ALL_STAGES)}")

        return stages

    def _load_prerequisites(self, stages):
        """
        Load data prerequisites when skipping early stages.

        When stages 0-1 are skipped, we still need sample_df and tRNA_data
        for later stages. Load them from existing files.
        """
        # Always need sample_df -- load from existing file or from config
        skip_early = not stages.intersection({'0a', '0b', '0c', '1'})
        if skip_early:
            # Try loading existing sample_df from output directory
            # Check both output root and data/ subdirectory
            sample_df_path = self.project_dir / 'sample_df.xlsx'
            if not sample_df_path.exists():
                sample_df_path = self.project_dir / 'data' / 'sample_df.xlsx'

            if sample_df_path.exists():
                self.log(f"Loading existing sample_df from {sample_df_path}")
                self.sample_df = pd.read_excel(sample_df_path)
                # Also load inp_file_df if it exists
                inp_file_path = self.project_dir / 'inp_file_df.xlsx'
                if not inp_file_path.exists():
                    inp_file_path = self.project_dir / 'data' / 'inp_file_df.xlsx'
                if inp_file_path.exists():
                    self.inp_file_df = pd.read_excel(inp_file_path)
                else:
                    self.inp_file_df = pd.DataFrame()
            else:
                # Fall back to loading from config
                self.log("No existing sample_df found, loading from config...")
                self.load_sample_info()

        # Stage 2 needs tRNA_data -- load from config if stage 1 was skipped
        if '2' in stages and '1' not in stages:
            if not hasattr(self, 'tRNA_data') or self.tRNA_data is None:
                self.log("Loading tRNA database for stats collection...")
                tRNA_database = self.config['tRNA_database']
                self.tRNA_data = read_tRNAdb_info(tRNA_database)

    def run(self, stages=None):
        """
        Run preprocessing pipeline.

        Args:
            stages: Set of stage identifiers to run, or None for all stages.
                    e.g. {'2', '3', '5'} to run only stats, charge, and QC.
        """
        if stages is None:
            stages = {'0a', '0b', '0c', '1', '2', '3', '4', '5', '6'}

        try:
            self.log("Starting tRNA-charge-seq preprocessing pipeline")
            self.log(f"Configuration: {self.config_file}")
            self.log(f"Project directory: {self.project_dir}")
            self.log(f"Parallel jobs: {self.n_jobs}")
            self.log(f"CPU budget: {self.n_jobs} jobs × {self.threads_per_job} threads = "
                     f"{self.n_jobs * self.threads_per_job} cores")
            self.log(f"Stages to run: {', '.join(sorted(stages))}")

            if self.sample_index is not None:
                self.log(f"Sample index: {self.sample_index} (single-sample mode)")
                if stages.intersection({'0a', '0b'}):
                    self.log("WARNING: --sample-index has no effect on stages 0a/0b "
                             "(they process all file pairs)", level="WARN")

            # Validate project directory for early stages
            if stages.intersection({'0a', '0b', '0c', '1'}):
                seq_dir = self.config.get('seq_dir', 'raw_fastq')
                raw_data_path = self.project_dir / 'data' / seq_dir
                if not raw_data_path.is_dir():
                    raise FileNotFoundError(
                        f"Raw data directory not found: {raw_data_path}\n"
                        f"--project-dir must point to the project root containing data/{seq_dir}/"
                    )

            # Setup directories
            self.setup_directories()

            # Load sample info for early stages
            if stages.intersection({'0a', '0b'}):
                # Stages 0a/0b: load fresh from config (sample list Excel)
                self.load_sample_info()
            elif stages.intersection({'0c', '1'}):
                # Stages 0c/1 without 0a/0b: need sample_df with N_total from 0b
                sample_df_path = self.project_dir / 'sample_df.xlsx'
                if sample_df_path.exists():
                    self.log(f"Loading existing sample_df from {sample_df_path}")
                    self.sample_df = pd.read_excel(sample_df_path)
                    inp_file_path = self.project_dir / 'inp_file_df.xlsx'
                    if inp_file_path.exists():
                        self.inp_file_df = pd.read_excel(inp_file_path)
                    else:
                        self.inp_file_df = pd.DataFrame()
                    self.log(f"Loaded {len(self.sample_df)} samples")
                    self.log(f"Input files: {len(self.inp_file_df)}")
                else:
                    self.log("WARNING: No sample_df.xlsx found — loading from config. "
                             "Stage 0c may fail if N_total column is missing.", level="WARN")
                    self.load_sample_info()

            # Load prerequisites for later stages when skipping early ones
            self._load_prerequisites(stages)

            # Run preprocessing stages (0-2)
            if '0a' in stages:
                t0 = time.time()
                self.stage_0a_adapter_removal()
                self.stage_timings['stage_0a'] = time.time() - t0
            if '0b' in stages:
                t0 = time.time()
                self.stage_0b_barcode_split()
                self.stage_timings['stage_0b'] = time.time() - t0
            if '0c' in stages:
                t0 = time.time()
                self.stage_0c_umi_trim()
                self.stage_timings['stage_0c'] = time.time() - t0
            if '1' in stages:
                t0 = time.time()
                self.stage_1_alignment()
                self.stage_timings['stage_1'] = time.time() - t0
            if '2' in stages:
                t0 = time.time()
                self.stage_2_stats_collection()
                self.stage_timings['stage_2'] = time.time() - t0

            # Run analysis stages (3-5)
            if '3' in stages:
                if self.config.get('run_charge_quantification', True):
                    t0 = time.time()
                    self.stage_3_charge_quantification()
                    self.stage_timings['stage_3'] = time.time() - t0
                else:
                    self.log("Skipping Stage 3: Charge Quantification (disabled in config)")

            if '4' in stages:
                if self.config.get('run_parquet_storage', False):
                    t0 = time.time()
                    self.stage_4_parquet_storage()
                    self.stage_timings['stage_4'] = time.time() - t0
                else:
                    self.log("Skipping Stage 4: Parquet Storage (disabled in config)")

            if '5' in stages:
                if self.config.get('run_qc_report', True):
                    t0 = time.time()
                    self.stage_5_qc_report()
                    self.stage_timings['stage_5'] = time.time() - t0
                else:
                    self.log("Skipping Stage 5: QC Report (disabled in config)")

            if '6' in stages:
                if self.config.get('run_modification_analysis', False):
                    t0 = time.time()
                    self.stage_6_modification_analysis()
                    self.stage_timings['stage_6'] = time.time() - t0
                else:
                    self.log("Skipping Stage 6: Modification Analysis (disabled in config)")

            # Collect file metrics and save
            if self.sample_index is not None:
                self._collect_file_metrics()
                self._save_metrics()
            else:
                self._collect_file_metrics()

            # Save outputs (skip in single-sample mode)
            if self.sample_index is None:
                self.save_outputs()
                self._save_metrics()
            else:
                self.log("Skipping save_outputs in single-sample mode")

            # Final summary
            self.status['end_time'] = datetime.now()
            duration = self.status['end_time'] - self.status['start_time']

            self.log("=" * 60)
            self.log("Pipeline Complete!")
            self.log("=" * 60)
            self.log(f"Duration: {duration}")
            self.log(f"Stages completed: {', '.join(self.status['stages_completed'])}")
            if self.stage_timings:
                self.log(f"Stage timings:")
                for stage, secs in sorted(self.stage_timings.items()):
                    m, s = divmod(int(secs), 60)
                    self.log(f"  {stage}: {m}m {s}s")
            self.log(f"Output files:")
            self.log(f"  - inp_file_df.xlsx")
            self.log(f"  - sample_df.xlsx")
            self.log(f"  - ALL_stats_aggregate.csv")

            # List charge results if generated
            if '3' in self.status['stages_completed']:
                self.log(f"  - charge_df_transcript.csv")
                self.log(f"  - charge_df_codon.csv")
                self.log(f"  - charge_df_aa.csv")
                self.log(f"  - charge_summary.csv")

            if '5' in self.status['stages_completed']:
                self.log(f"  - QC_summary.csv")
                self.log(f"  - QC_report.html")

            self.log(f"  - computing_metrics.csv")
            self.log("")
            if '3' not in self.status['stages_completed']:
                self.log("Next step: Run charge quantification (Stage 3)")
                self.log("  Re-run with --stages 3 or:")
                self.log("  python -m trnaseq quantify -i ALL_stats_aggregate.csv -o charge.csv")

        except Exception as e:
            self.log(f"Pipeline failed: {str(e)}", level="ERROR")
            raise


def main():
    parser = argparse.ArgumentParser(
        description="Unified tRNA-charge-seq preprocessing pipeline (Stages 0-5)"
    )
    parser.add_argument(
        '--config', required=True,
        help='YAML configuration file'
    )
    parser.add_argument(
        '--project-dir', required=False, default=None,
        help='Project directory (must contain data/raw_fastq/ for stages 0-1)'
    )
    parser.add_argument(
        '--output-dir', required=False, default=None, dest='output_dir_deprecated',
        help='Deprecated: use --project-dir instead'
    )
    parser.add_argument(
        '--n-jobs', type=int, default=4,
        help='Number of parallel jobs (default: 4)'
    )
    parser.add_argument(
        '--skip-charge', action='store_true',
        help='Skip charge quantification (Stage 3)'
    )
    parser.add_argument(
        '--parquet', action='store_true',
        help='Save to Parquet format (Stage 4)'
    )
    parser.add_argument(
        '--stages', type=str, default='all',
        help='Stages to run: "all", "0-2", "3-5", or specific like "0a,0b,1,3"'
    )
    parser.add_argument(
        '--sample-index', type=int, default=None,
        help='Process only sample at this 0-based index (for SLURM array jobs)'
    )
    parser.add_argument(
        '--threads-per-job', type=int, default=None,
        help='Threads per subprocess (AdapterRemoval/SWIPE). '
             'Default: config threads_per_job or 2'
    )

    args = parser.parse_args()

    # Resolve --project-dir / --output-dir (deprecated)
    project_dir = args.project_dir or args.output_dir_deprecated
    if project_dir is None:
        parser.error('--project-dir is required')
    if args.output_dir_deprecated and not args.project_dir:
        print('Warning: --output-dir is deprecated, use --project-dir instead',
              file=sys.stderr)

    # Run pipeline
    pipeline = PreprocessingPipeline(
        config_file=args.config,
        project_dir=project_dir,
        n_jobs=args.n_jobs,
        sample_index=args.sample_index,
        threads_per_job=args.threads_per_job,
    )

    # Override config with CLI arguments
    if args.skip_charge:
        pipeline.config['run_charge_quantification'] = False
    if args.parquet:
        pipeline.config['run_parquet_storage'] = True

    # Parse stages
    stages = PreprocessingPipeline.parse_stages(args.stages)

    pipeline.run(stages=stages)


if __name__ == '__main__':
    main()
