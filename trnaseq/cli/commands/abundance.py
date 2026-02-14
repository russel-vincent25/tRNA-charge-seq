"""
CLI command: trnaseq abundance
================================

Run differential tRNA abundance analysis using pyDESeq2.

Usage:
    python -m trnaseq abundance \
        -i ALL_stats_aggregate.csv \
        --sample-df sample_df.xlsx \
        --level aa \
        --control WT \
        -o results/
"""

import argparse
import sys


def add_abundance_parser(subparsers):
    """Register the 'abundance' subcommand."""
    parser = subparsers.add_parser(
        'abundance',
        help='Differential tRNA abundance analysis (pyDESeq2)'
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help='Path to ALL_stats_aggregate.csv'
    )
    parser.add_argument(
        '--sample-df', required=True,
        help='Path to sample_df.xlsx (must have sample_name_unique and sample_name columns)'
    )
    parser.add_argument(
        '--level', choices=['transcript', 'codon', 'aa'], default='aa',
        help='Aggregation level (default: aa)'
    )
    parser.add_argument(
        '--control', default=None,
        help='Control group name (must match a sample_name value). '
             'Defaults to first group alphabetically.'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Output directory for results'
    )
    parser.set_defaults(func=run_abundance)


def run_abundance(args):
    """Execute the abundance command."""
    try:
        from trnaseq.abundance import DifferentialAbundance
    except ImportError:
        print("ERROR: pyDESeq2 is required for differential abundance analysis.",
              file=sys.stderr)
        print("Install with: pip install pydeseq2", file=sys.stderr)
        sys.exit(1)

    print(f"Input:   {args.input}")
    print(f"Samples: {args.sample_df}")
    print(f"Level:   {args.level}")
    print(f"Control: {args.control or '(auto)'}")
    print(f"Output:  {args.output}")
    print()

    da = DifferentialAbundance(
        stats_csv=args.input,
        sample_df=args.sample_df,
        level=args.level,
        control_group=args.control,
    )

    print(f"Count matrix: {da.count_matrix.shape[0]} samples x {da.count_matrix.shape[1]} features")
    print(f"Control group: {da.control_group}")
    print()

    print("Running DESeq2...")
    results = da.run_deseq2()
    print(f"Results: {len(results)} comparisons")

    if not results.empty:
        sig = results['padj'] < 0.05
        print(f"Significant (padj < 0.05): {sig.sum()}")

    print()
    print("Exporting results...")
    da.export_results(args.output)
    print(f"Done! Results saved to {args.output}/")
