"""
CLI command for building tRNA reference databases from MODOMICS.

Usage:
    python -m trnaseq build-reference --organism human --output-dir tRNA_database/human/
    python -m trnaseq build-reference --organism ecoli --output-dir tRNA_database/ecoli_modomics/ --full-table
"""

import argparse
import sys


def add_parser(subparsers):
    """Add the build-reference subcommand to the CLI."""
    parser = subparsers.add_parser(
        'build-reference',
        help='Build tRNA reference FASTA from MODOMICS isodecoder sequences',
        description=(
            'Generate a non-redundant tRNA reference FASTA and modification '
            'table from MODOMICS. Each entry is a unique isodecoder (mature '
            'tRNA transcript), eliminating gene copy redundancy that causes '
            'multi-mapping in eukaryotic genomes.'
        ),
    )
    parser.add_argument(
        '--organism', '-g', required=True,
        help='Organism name (e.g. ecoli, human, mouse, "Saccharomyces cerevisiae")',
    )
    parser.add_argument(
        '--output-dir', '-o', required=True,
        help='Output directory for FASTA and modification CSV files',
    )
    parser.add_argument(
        '--name-prefix', '-p', default=None,
        help='Prefix for FASTA sequence names (default: organism name)',
    )
    parser.add_argument(
        '--full-table', action='store_true',
        help='Also generate a complete per-position modification table '
             '(like all_tRNA_pos_mod_info.csv)',
    )
    parser.add_argument(
        '--no-api', action='store_true',
        help='Skip MODOMICS API and use CSV fallback only',
    )
    parser.set_defaults(func=run)


def run(args):
    """Execute the build-reference command."""
    from trnaseq.modifications.reference_builder import (
        build_modomics_reference,
        build_full_position_table,
    )

    print(f'Building MODOMICS reference for: {args.organism}')
    print(f'Output directory: {args.output_dir}')

    try:
        fasta_path, mod_csv_path = build_modomics_reference(
            organism=args.organism,
            output_dir=args.output_dir,
            use_api=not args.no_api,
            name_prefix=args.name_prefix,
        )
        print(f'  FASTA:  {fasta_path}')
        print(f'  Mods:   {mod_csv_path}')

        # Count sequences
        n_seqs = sum(1 for line in open(fasta_path) if line.startswith('>'))
        n_mods = sum(1 for _ in open(mod_csv_path)) - 1  # subtract header
        print(f'  {n_seqs} isodecoder sequences, {n_mods} modification sites')

        if args.full_table:
            safe_org = args.organism.lower().replace(' ', '_')
            table_path = build_full_position_table(
                organism=args.organism,
                output_path=f'{args.output_dir}/{safe_org}_all_pos_mod_info.csv',
                use_api=not args.no_api,
                name_prefix=args.name_prefix,
            )
            n_rows = sum(1 for _ in open(table_path)) - 1
            print(f'  Table:  {table_path} ({n_rows} positions)')

    except ValueError as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f'Failed: {e}', file=sys.stderr)
        sys.exit(1)

    print('Done.')
