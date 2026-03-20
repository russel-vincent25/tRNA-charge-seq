"""
CLI command: trnaseq fragments
================================

Per-tRNA fragment classification, RT drop-off profiling, and fragment
length distributions.

Usage:
    python -m trnaseq fragments \
        --stats-dir data/stats_collection/ \
        --min-reads 10 \
        --csv \
        -o results/fragments/
"""

import argparse


def add_fragments_parser(subparsers):
    """Register the 'fragments' subcommand."""
    parser = subparsers.add_parser(
        'fragments',
        help='Fragment classification and RT drop-off profiling'
    )
    parser.add_argument('--stats-dir', required=True,
                        help='Directory containing *_stats.csv.bz2 files')
    parser.add_argument('--sample-df', default=None,
                        help='Path to sample_df.xlsx (optional; auto-detects samples if omitted)')
    parser.add_argument('--output-dir', '-o', required=True,
                        help='Output directory for results')
    parser.add_argument('--min-reads', type=int, default=10,
                        help='Minimum reads per tRNA to include (default: 10)')
    parser.add_argument('--csv', action='store_true',
                        help='Also write CSV copies (default: parquet only)')
    parser.set_defaults(func=run_fragments)


def run_fragments(args):
    """Execute the fragments command."""
    from pathlib import Path
    import pandas as pd
    from trnaseq.fragments import FragmentAnalyser

    stats_dir = Path(args.stats_dir)
    output_dir = Path(args.output_dir)

    # Determine sample names
    sample_names = None
    if args.sample_df:
        sample_df = pd.read_excel(args.sample_df)
        if 'sample_name_unique' in sample_df.columns:
            sample_names = sample_df['sample_name_unique'].tolist()
        else:
            sample_names = sample_df.iloc[:, 0].tolist()

    print(f"Stats dir:   {stats_dir}")
    print(f"Output:      {output_dir}")
    print(f"Min reads:   {args.min_reads}")
    print()

    analyser = FragmentAnalyser(
        stats_dir=stats_dir,
        sample_names=sample_names,
        min_reads=args.min_reads,
    )

    print(f"Found {len(analyser.sample_names)} samples")
    print()

    analyser.run()
    analyser.export(output_dir, write_csv=args.csv)

    # Print summary
    summary = analyser._summary
    if summary is not None and not summary.empty:
        print(f"\nFragment summary ({len(summary)} samples):")
        for _, row in summary.iterrows():
            print(f"  {row['sample_name_unique']}: "
                  f"{row['pct_full_length']:.1f}% full-length, "
                  f"{row['pct_rt_dropoff']:.1f}% RT drop-off, "
                  f"{row['pct_5p_tRF']:.1f}% 5'tRF, "
                  f"{row['pct_degraded']:.1f}% degraded")

    print(f"\nDone! Results saved to {output_dir}/")
