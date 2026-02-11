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
- Stage 5: Alignment QC Reports (optional)

Output:
- inp_file_df.xlsx (input file summary)
- sample_df.xlsx (sample information with QC metrics)
- ALL_stats_aggregate.csv (combined statistics for charge quantification)
- charge_analysis/ (charge quantification results, if enabled)
- parquet_data/ (parquet-format data, if enabled)
- qc_reports/ (alignment QC reports, if enabled)

Usage:
    python scripts/run_preprocessing.py \
        --config config.yaml \
        --output-dir output/ \
        --n-jobs 8 \
        --parquet \
        --qc-reports

Author: Based on projects/example/process_data.ipynb
Date: 2026-02-11
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
    from trnaseq.visualization.alignment_viewer import AlignmentViewer
    VIEWER_AVAILABLE = True
except ImportError:
    VIEWER_AVAILABLE = False


class PreprocessingPipeline:
    """
    Unified preprocessing pipeline for tRNA-charge-seq

    Runs stages 0a-5 and outputs key files for downstream analysis.
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
            'charge_dir': 'charge_analysis',
            'parquet_dir': 'parquet_data',
            'qc_dir': 'qc_reports',
        }

        # Create data directory
        (self.output_dir / 'data').mkdir(exist_ok=True)

        # Create subdirectories under data/
        for dir_name in ['AdapterRemoval', 'BC_split', 'UMI_trimmed',
                         'SWalign', 'stats_collection']:
            (self.output_dir / 'data' / dir_name).mkdir(exist_ok=True)

        # Create analysis directories
        for dir_name in ['charge_analysis', 'parquet_data', 'qc_reports']:
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
        stats_file = self.output_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
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
            charge_dir = self.output_dir / 'charge_analysis'
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
            parquet_dir = self.output_dir / 'parquet_data'
            store = tRNAseqDataStore(str(parquet_dir))

            compression = self.config.get('parquet_compression', 'snappy')
            self.log(f"  Compression: {compression}")

            # Convert stats CSV to Parquet
            stats_file = self.output_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
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
            charge_dir = self.output_dir / 'charge_analysis'
            if charge_dir.exists():
                for csv_file in charge_dir.glob('charge_df_*.csv'):
                    self.log(f"  Converting {csv_file.name} to Parquet...")
                    df = pd.read_csv(csv_file)
                    parquet_file = store.save_charge_data(
                        df,
                        compression=compression
                    )
                    # Rename to match CSV name
                    new_name = parquet_dir / 'charge' / csv_file.with_suffix('.parquet').name
                    parquet_file.rename(new_name)
                    self.log(f"    Saved: {new_name.name}")

            self.status['stages_completed'].append('4')
            self.log("  Parquet conversion complete!")

        except Exception as e:
            self.log(f"ERROR in Parquet storage: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def stage_5_alignment_qc(self):
        """Stage 5: Alignment QC Reports (optional)"""
        self.log("=" * 60)
        self.log("Stage 5: Alignment QC Reports")
        self.log("=" * 60)

        if not VIEWER_AVAILABLE:
            self.log("WARNING: AlignmentViewer not available. Skipping stage 5.", level="WARN")
            self.log("  Install with: pip install -e . (from repo root)", level="WARN")
            return

        try:
            qc_dir = self.output_dir / 'qc_reports'
            qc_dir.mkdir(exist_ok=True)

            # Get QC settings
            qc_format = self.config.get('qc_report_format', 'html')
            top_n = self.config.get('qc_top_trnas', 20)

            # Find alignment JSON files
            align_dir = self.output_dir / 'data' / 'SWalign'
            json_files = list(align_dir.glob('*_SWalign.json.bz2'))

            if not json_files:
                self.log("  No alignment JSON files found", level="WARN")
                return

            self.log(f"  Found {len(json_files)} alignment files")
            self.log(f"  Generating QC reports for top {top_n} tRNAs per sample")

            for json_file in json_files:
                sample_name = json_file.name.replace('_SWalign.json.bz2', '')
                self.log(f"  Processing {sample_name}...")

                try:
                    viewer = AlignmentViewer(str(json_file))

                    # List top tRNAs by read count
                    trna_df = viewer.list_trnas(min_reads=10)

                    if trna_df.empty:
                        self.log(f"    No tRNAs with >= 10 reads", level="WARN")
                        continue

                    # Save tRNA list
                    trna_list_file = qc_dir / f'{sample_name}_trna_list.csv'
                    trna_df.to_csv(trna_list_file, index=False)
                    self.log(f"    Saved tRNA list: {trna_list_file.name}")

                    # Generate reports for top N tRNAs
                    top_trnas = trna_df.head(top_n)['tRNA'].tolist()

                    for i, trna_id in enumerate(top_trnas[:5], 1):  # Limit to 5 for speed
                        self.log(f"    Generating report for {trna_id} ({i}/5)...")

                        if qc_format == 'html':
                            output_file = qc_dir / f'{sample_name}_{trna_id}_report.html'
                            viewer.create_html_report(trna_id, output=str(output_file))
                        else:
                            output_file = qc_dir / f'{sample_name}_{trna_id}_coverage.png'
                            viewer.plot_coverage(trna_id, output=str(output_file))

                except Exception as e:
                    self.log(f"    ERROR processing {sample_name}: {str(e)}", level="ERROR")
                    continue

            self.status['stages_completed'].append('5')
            self.log("  QC report generation complete!")

        except Exception as e:
            self.log(f"ERROR in QC reports: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

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

        # Copy charge analysis results to output root
        if hasattr(self, 'charge_results') and self.charge_results:
            charge_dir = self.output_dir / 'charge_analysis'
            for level, df in self.charge_results.items():
                output_file = self.output_dir / f'charge_df_{level}.csv'
                import shutil
                shutil.copy(
                    charge_dir / f'charge_df_{level}.csv',
                    output_file
                )
                self.log(f"  Saved: {output_file.name}")

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

            # Run preprocessing stages (0-2)
            self.stage_0a_adapter_removal()
            self.stage_0b_barcode_split()
            self.stage_0c_umi_trim()
            self.stage_1_alignment()
            self.stage_2_stats_collection()

            # Run optional analysis stages (3-5)
            if self.config.get('run_charge_quantification', True):
                self.stage_3_charge_quantification()
            else:
                self.log("Skipping Stage 3: Charge Quantification (disabled in config)")

            if self.config.get('run_parquet_storage', False):
                self.stage_4_parquet_storage()
            else:
                self.log("Skipping Stage 4: Parquet Storage (disabled in config)")

            if self.config.get('run_alignment_qc', False):
                self.stage_5_alignment_qc()
            else:
                self.log("Skipping Stage 5: Alignment QC (disabled in config)")

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

            # List charge results if generated
            if '3' in self.status['stages_completed']:
                self.log(f"  - charge_df_transcript.csv")
                self.log(f"  - charge_df_codon.csv")
                self.log(f"  - charge_df_aa.csv")
                self.log(f"  - charge_summary.csv")

            self.log("")
            if '3' not in self.status['stages_completed']:
                self.log("Next step: Run charge quantification (Stage 3)")
                self.log("  Re-run with charge quantification enabled in config, or:")
                self.log("  python -m trnaseq.cli.commands.quantify ALL_stats_aggregate.csv")

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
        '--output-dir', required=True,
        help='Output directory'
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
        '--qc-reports', action='store_true',
        help='Generate alignment QC reports (Stage 5)'
    )
    parser.add_argument(
        '--stages', type=str, default='all',
        help='Stages to run: "all", "0-2", "3-5", or specific like "0a,0b,1,3"'
    )

    args = parser.parse_args()

    # Run pipeline
    pipeline = PreprocessingPipeline(
        config_file=args.config,
        output_dir=args.output_dir,
        n_jobs=args.n_jobs
    )

    # Override config with CLI arguments
    if args.skip_charge:
        pipeline.config['run_charge_quantification'] = False
    if args.parquet:
        pipeline.config['run_parquet_storage'] = True
    if args.qc_reports:
        pipeline.config['run_alignment_qc'] = True

    # Note: --stages argument is parsed but not yet implemented
    # This would require refactoring run() to support selective stage execution
    if args.stages != 'all':
        print(f"WARNING: --stages option not yet implemented. Running all stages.")
        print(f"  Requested stages: {args.stages}")

    pipeline.run()


if __name__ == '__main__':
    main()
