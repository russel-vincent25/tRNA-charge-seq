"""
Alignment View CLI Command

Wraps AlignmentViewer for on-demand alignment visualization.
Replaces the old per-tRNA QC reports that ran inside the pipeline.

Usage:
    python -m trnaseq view --json data/SWalign/sample_SWalign.json.bz2 --trna tRNA-Ala-TGC-1
    python -m trnaseq view --json data/SWalign/sample_SWalign.json.bz2 --list
"""

import argparse
import sys
from pathlib import Path


def add_view_parser(subparsers):
    """Add the 'view' subcommand to the argument parser."""
    parser = subparsers.add_parser(
        'view',
        help='View tRNA alignment coverage plots and reports',
        description='On-demand alignment visualization from SWIPE JSON files'
    )

    parser.add_argument(
        '--json', required=True, type=str,
        help='Path to {sample}_SWalign.json.bz2 file'
    )

    # Mutually exclusive: list tRNAs or view one
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--trna', type=str, default=None,
        help='tRNA ID to visualize (e.g. tRNA-Ala-TGC-1)'
    )
    group.add_argument(
        '--list', action='store_true', dest='list_trnas',
        help='List all tRNAs in the alignment file with read counts'
    )

    parser.add_argument(
        '--output', '-o', type=str, default=None,
        help='Output file path (.png or .html). Default: auto-named in current dir'
    )
    parser.add_argument(
        '--format', choices=['html', 'png'], default='html',
        help='Report format (default: html)'
    )
    parser.add_argument(
        '--min-reads', type=int, default=10,
        help='Minimum reads to include in --list output (default: 10)'
    )

    parser.set_defaults(func=run_view)


def run_view(args):
    """Execute the view command."""
    json_path = Path(args.json)
    if not json_path.exists():
        print(f'Error: JSON file not found: {args.json}', file=sys.stderr)
        sys.exit(1)

    try:
        from trnaseq.visualization.alignment_viewer import AlignmentViewer
    except ImportError:
        print('Error: AlignmentViewer not available. Check trnaseq installation.',
              file=sys.stderr)
        sys.exit(1)

    viewer = AlignmentViewer(str(json_path))

    if args.list_trnas:
        # List mode
        trna_df = viewer.list_trnas(min_reads=args.min_reads)
        if trna_df.empty:
            print(f'No tRNAs with >= {args.min_reads} reads')
        else:
            print(trna_df.to_string(index=False))

        # Optionally save to CSV
        if args.output:
            trna_df.to_csv(args.output, index=False)
            print(f'\nSaved to {args.output}')
        return

    # Single tRNA view mode
    trna_id = args.trna
    sample_name = json_path.stem.replace('_SWalign.json', '').replace('.bz2', '')

    if args.output:
        output = args.output
    elif args.format == 'html':
        output = f'{sample_name}_{trna_id}_report.html'
    else:
        output = f'{sample_name}_{trna_id}_coverage.png'

    if args.format == 'html':
        viewer.create_html_report(trna_id, output=output)
    else:
        viewer.plot_coverage(trna_id, output=output)

    print(f'Output: {output}')
