#!/usr/bin/env python3
"""
Unified Preprocessing Pipeline (Stages 0-2)
============================================

Combines all preprocessing steps into a single script:
- Stage 0a: Adapter Removal & Merging
- Stage 0b: Barcode Splitting
- Stage 0c: UMI Trimming
- Stage 1: SWIPE Alignment
- Stage 2: Stats Collection

Output:
- inp_file_df.xlsx (input file summary)
- sample_df.xlsx (sample information with QC metrics)
- ALL_stats_aggregate.csv (combined statistics for charge quantification)

Usage:
    python scripts/run_preprocessing.py \
        --config config.yaml \
        --output-dir output/ \
        --n-jobs 8

Author: Based on projects/example/process_data.ipynb
Date: 2026-02-10
"""

import os
import sys
from pathlib import Path
from datetime import datetime
import argparse
import yaml
import pandas as pd
import numpy as np

# Add repo to path
repo_path = Path(__file__).parent.parent
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


class PreprocessingPipeline:
    """
    Unified preprocessing pipeline for tRNA-charge-seq

    Runs stages 0a-2 and outputs key files for downstream analysis.
    """

    def __init__(self, config_file, output_dir, n_jobs=4):
        """
        Initialize pipeline

        Parameters:
            config_file: YAML configuration file
            output_dir: Output directory for results
            n_jobs: Number of parallel jobs
        """
        self.config_file = Path(config_file)
        self.output_dir = Path(output_dir)
        self.n_jobs = n_jobs
        self.log_file = self.output_dir / "preprocessing.log"

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load configuration
        with open(self.config_file, 'r') as f:
            self.config = yaml.safe_load(f)

        # Initialize status tracking
        self.status = {
            'start_time': datetime.now(),
            'stages_completed': []
        }

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
            'NBdir': str(self.output_dir),
            'data_dir': 'data',
            'seq_dir': self.config.get('seq_dir', 'raw_fastq'),
            'AdapterRemoval_dir': 'AdapterRemoval',
            'BC_dir': 'BC_split',
            'UMI_dir': 'UMI_trimmed',
            'align_dir': 'SWalign',
            'stats_dir': 'stats_collection',
        }

        # Create subdirectories
        for dir_name in ['data', 'AdapterRemoval', 'BC_split', 'UMI_trimmed',
                         'SWalign', 'stats_collection']:
            (self.output_dir / dir_name).mkdir(exist_ok=True)

    def stage_0a_adapter_removal(self):
        """Stage 0a: Adapter Removal & Merging"""
        self.log("=" * 60)
        self.log("Stage 0a: Adapter Removal & Merging")
        self.log("=" * 60)

        MIN_READ_LEN = self.config.get('min_read_len', 39)

        # Run AdapterRemoval
        AR_obj = AR_merge(
            self.dir_dict, self.inp_file_df, MIN_READ_LEN,
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

        UMItrim_obj = UMI_trim(
            self.dir_dict, self.sample_df,
            overwrite_dir=self.config.get('overwrite', True),
            downsample_percentile=downsample_percentile,
            downsample_absolute=downsample_absolute
        )
        self.sample_df = UMItrim_obj.run_parallel(n_jobs=self.n_jobs)

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
        align_obj = SWIPE_align(
            self.dir_dict, tRNA_database, self.sample_df,
            SWIPE_score_mat,
            gap_penalty=self.config.get('gap_penalty', 6),
            extension_penalty=self.config.get('extension_penalty', 3),
            min_score_align=MIN_SCORE_ALIGN,
            common_seqs=common_seqs,
            overwrite_dir=self.config.get('overwrite', True)
        )
        self.sample_df = align_obj.run_parallel(n_jobs=self.n_jobs)

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

        # Run stats collection
        stats_obj = STATS_collection(
            self.dir_dict, self.tRNA_data, self.sample_df,
            common_seqs=common_seqs,
            overwrite_dir=self.config.get('overwrite', True)
        )
        self.stats_df = stats_obj.run_parallel(n_jobs=self.n_jobs)

        self.log(f"  Collected stats for {len(self.stats_df)} entries")

        self.status['stages_completed'].append('2')

    def save_outputs(self):
        """Save output files"""
        self.log("=" * 60)
        self.log("Saving outputs")
        self.log("=" * 60)

        # Save inp_file_df
        inp_file_output = self.output_dir / 'inp_file_df.xlsx'
        self.inp_file_df.to_excel(inp_file_output, index=False)
        self.log(f"  Saved: {inp_file_output}")

        # Save sample_df
        sample_df_output = self.output_dir / 'sample_df.xlsx'
        self.sample_df.to_excel(sample_df_output, index=False)
        self.log(f"  Saved: {sample_df_output}")

        # ALL_stats_aggregate.csv should already be in stats_collection/
        stats_file = self.output_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
        if stats_file.exists():
            # Copy to output root for easy access
            import shutil
            shutil.copy(stats_file, self.output_dir / 'ALL_stats_aggregate.csv')
            self.log(f"  Saved: ALL_stats_aggregate.csv")
        else:
            self.log(f"  WARNING: ALL_stats_aggregate.csv not found", level="WARN")

    def run(self):
        """Run complete preprocessing pipeline"""
        try:
            self.log("Starting tRNA-charge-seq preprocessing pipeline")
            self.log(f"Configuration: {self.config_file}")
            self.log(f"Output directory: {self.output_dir}")
            self.log(f"Parallel jobs: {self.n_jobs}")

            # Setup
            self.load_sample_info()
            self.setup_directories()

            # Run stages
            self.stage_0a_adapter_removal()
            self.stage_0b_barcode_split()
            self.stage_0c_umi_trim()
            self.stage_1_alignment()
            self.stage_2_stats_collection()

            # Save outputs
            self.save_outputs()

            # Final summary
            self.status['end_time'] = datetime.now()
            duration = self.status['end_time'] - self.status['start_time']

            self.log("=" * 60)
            self.log("Pipeline Complete!")
            self.log("=" * 60)
            self.log(f"Duration: {duration}")
            self.log(f"Stages completed: {', '.join(self.status['stages_completed'])}")
            self.log(f"Output files:")
            self.log(f"  - inp_file_df.xlsx")
            self.log(f"  - sample_df.xlsx")
            self.log(f"  - ALL_stats_aggregate.csv")
            self.log("")
            self.log("Next step: Run charge quantification (Stage 3)")
            self.log("  python -m trnaseq.cli.commands.quantify ALL_stats_aggregate.csv")

        except Exception as e:
            self.log(f"Pipeline failed: {str(e)}", level="ERROR")
            raise


def main():
    parser = argparse.ArgumentParser(
        description="Unified tRNA-charge-seq preprocessing pipeline (Stages 0-2)"
    )
    parser.add_argument(
        '--config', required=True,
        help='YAML configuration file'
    )
    parser.add_argument(
        '--output-dir', required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--n-jobs', type=int, default=4,
        help='Number of parallel jobs (default: 4)'
    )

    args = parser.parse_args()

    # Run pipeline
    pipeline = PreprocessingPipeline(
        config_file=args.config,
        output_dir=args.output_dir,
        n_jobs=args.n_jobs
    )
    pipeline.run()


if __name__ == '__main__':
    main()
