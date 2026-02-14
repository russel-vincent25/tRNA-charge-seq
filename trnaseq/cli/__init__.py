"""
tRNA-seq Command Line Interface

Usage:
    python -m trnaseq pipeline --config config.yaml --project-dir /path/to/project
    python -m trnaseq quantify -i stats.csv -o charge.csv
    python -m trnaseq view --json sample_SWalign.json.bz2 --trna tRNA-Ala-TGC-1
"""

import argparse
import sys


def main(argv=None):
    """Main CLI dispatcher for trnaseq subcommands."""
    parser = argparse.ArgumentParser(
        prog='trnaseq',
        description='tRNA-seq analysis toolkit'
    )
    subparsers = parser.add_subparsers(dest='command')

    # --- pipeline subcommand ---
    pipeline_parser = subparsers.add_parser(
        'pipeline',
        help='Run preprocessing pipeline (stages 0-5)'
    )
    pipeline_parser.add_argument('--config', required=True, help='YAML configuration file')
    pipeline_parser.add_argument(
        '--project-dir', required=False, default=None,
        help='Project directory (must contain data/raw_fastq/ for stages 0-1)'
    )
    pipeline_parser.add_argument(
        '--output-dir', required=False, default=None, dest='output_dir_deprecated',
        help='Deprecated: use --project-dir instead'
    )
    pipeline_parser.add_argument('--n-jobs', type=int, default=4, help='Parallel jobs (default: 4)')
    pipeline_parser.add_argument('--skip-charge', action='store_true', help='Skip charge quantification')
    pipeline_parser.add_argument('--parquet', action='store_true', help='Save to Parquet format')
    pipeline_parser.add_argument(
        '--stages', type=str, default='all',
        help='Stages to run: "all", "0-2", "3-5", or "0a,0b,1,3"'
    )
    pipeline_parser.add_argument(
        '--sample-index', type=int, default=None,
        help='Process only sample at this 0-based index (for SLURM array jobs)'
    )

    # --- quantify subcommand ---
    from trnaseq.cli.commands.quantify import add_quantify_parser
    add_quantify_parser(subparsers)

    # --- view subcommand ---
    from trnaseq.cli.commands.view import add_view_parser
    add_view_parser(subparsers)

    # --- abundance subcommand ---
    from trnaseq.cli.commands.abundance import add_abundance_parser
    add_abundance_parser(subparsers)

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == 'pipeline':
        _run_pipeline(args, pipeline_parser)
    elif args.command in ('quantify', 'view', 'abundance'):
        args.func(args)


def _run_pipeline(args, parser):
    """Handle the pipeline subcommand."""
    from trnaseq.cli.run_pipeline import PreprocessingPipeline

    project_dir = args.project_dir or args.output_dir_deprecated
    if project_dir is None:
        parser.error('--project-dir is required')
    if args.output_dir_deprecated and not args.project_dir:
        print('Warning: --output-dir is deprecated, use --project-dir instead',
              file=sys.stderr)

    pipeline = PreprocessingPipeline(
        config_file=args.config,
        project_dir=project_dir,
        n_jobs=args.n_jobs,
        sample_index=getattr(args, 'sample_index', None),
    )

    if args.skip_charge:
        pipeline.config['run_charge_quantification'] = False
    if args.parquet:
        pipeline.config['run_parquet_storage'] = True

    stages = PreprocessingPipeline.parse_stages(args.stages)
    pipeline.run(stages=stages)
