"""
Charge Quantification CLI Command

This module provides the 'quantify' command for calculating tRNA charge levels
from alignment statistics CSV files.

Usage:
    trnaseq quantify --input ALL_stats_aggregate.csv --output charge_data.csv
"""

import argparse
import sys
from pathlib import Path
from trnaseq.charge import ChargeQuantifier


def add_quantify_parser(subparsers):
    """
    Add the 'quantify' subcommand to the argument parser.

    Args:
        subparsers: ArgumentParser subparsers object
    """
    parser = subparsers.add_parser(
        'quantify',
        help='Quantify tRNA charge levels from alignment statistics',
        description='Calculate tRNA charging percentages from ALL_stats_aggregate.csv files'
    )

    # Required arguments
    parser.add_argument(
        '-i', '--input',
        required=True,
        type=str,
        help='Path to ALL_stats_aggregate.csv file'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        type=str,
        help='Path to output CSV file for charge data'
    )

    # Optional arguments
    parser.add_argument(
        '-l', '--level',
        type=str,
        choices=['transcript', 'codon', 'aa'],
        default='transcript',
        help='Annotation level for charge quantification (default: transcript)'
    )

    parser.add_argument(
        '--charge-count',
        type=str,
        choices=['count', 'UMIcount'],
        default='count',
        help='Column to use for charge calculation (default: count)'
    )

    parser.add_argument(
        '--rpm-count',
        type=str,
        choices=['count', 'UMIcount'],
        default='UMIcount',
        help='Column to use for RPM calculation (default: UMIcount)'
    )

    parser.add_argument(
        '--exclude-gaps',
        action='store_true',
        help='Exclude alignments with gaps'
    )

    parser.add_argument(
        '--exclude-low-score',
        action='store_true',
        help='Exclude alignments with fmax_score < 0.9'
    )

    parser.add_argument(
        '--include-synthetic',
        action='store_true',
        help='Include synthetic control tRNAs in output'
    )

    parser.add_argument(
        '--exclude-mito',
        action='store_true',
        help='Exclude mitochondrial tRNAs from output'
    )

    parser.add_argument(
        '--summary',
        type=str,
        default=None,
        help='Optional path to output summary statistics CSV'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )

    parser.set_defaults(func=run_quantify)


def run_quantify(args):
    """
    Execute the charge quantification command.

    Args:
        args: Parsed command-line arguments
    """
    # Validate input file exists
    input_path = Path(args.input)
    if not input_path.exists():
        print(f'Error: Input file not found: {args.input}', file=sys.stderr)
        sys.exit(1)

    # Print header
    if args.verbose:
        print('=' * 60)
        print('tRNA Charge Quantification')
        print('=' * 60)
        print(f'Input file: {args.input}')
        print(f'Output file: {args.output}')
        print(f'Annotation level: {args.level}')
        print(f'Charge count column: {args.charge_count}')
        print(f'RPM count column: {args.rpm_count}')
        print('-' * 60)

    try:
        # Initialize ChargeQuantifier
        if args.verbose:
            print('Loading alignment statistics...')

        quantifier = ChargeQuantifier(
            stats_csv=args.input,
            charge_count=args.charge_count,
            RPM_count=args.rpm_count,
            excl_align_gap=args.exclude_gaps,
            excl_09_fmax=args.exclude_low_score
        )

        if args.verbose:
            print(f'Loaded {len(quantifier.stats_df)} rows from input file')
            print('Calculating charge levels...')

        # Quantify charge
        charge_data = quantifier.quantify_all(
            level=args.level,
            include_synthetic=args.include_synthetic,
            include_mito=not args.exclude_mito
        )

        if args.verbose:
            print(f'Calculated charge for {len(charge_data)} entries')

        # Export to CSV
        if args.verbose:
            print(f'Writing charge data to {args.output}...')

        charge_data.to_csv(args.output, index=False)

        # Generate summary statistics if requested
        if args.summary:
            if args.verbose:
                print(f'Generating summary statistics...')

            summary = quantifier.get_summary_statistics(level=args.level)
            summary.to_csv(args.summary, index=False)

            if args.verbose:
                print(f'Wrote summary statistics to {args.summary}')

        # Print success message
        if args.verbose:
            print('-' * 60)
            print('Charge quantification completed successfully!')
            print('=' * 60)
        else:
            print(f'Charge data written to {args.output}')

    except Exception as e:
        print(f'Error during charge quantification: {str(e)}', file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def main():
    """
    Standalone entry point: python -m trnaseq.cli.commands.quantify
    """
    parser = argparse.ArgumentParser(
        description='Quantify tRNA charge levels from alignment statistics'
    )
    subparsers = parser.add_subparsers(dest='command')
    add_quantify_parser(subparsers)
    args = parser.parse_args()

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
