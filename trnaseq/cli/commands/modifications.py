"""
CLI command: trnaseq modifications
====================================

Run per-position RT signature analysis and modification calling.

Usage:
    python -m trnaseq modifications \
        --json-dir data/SWalign/ \
        --reference tRNA_database/human/human-tRNAs.fa \
        --output-dir modification_analysis/ \
        --organism "Homo sapiens" \
        --discover-novel
"""

import argparse
import sys


def add_modifications_parser(subparsers):
    """Register the 'modifications' subcommand."""
    parser = subparsers.add_parser(
        'modifications',
        help='Per-position RT analysis and modification calling'
    )
    parser.add_argument('--json-dir', required=True,
                        help='Directory containing *_SWalign.json.bz2 files')
    parser.add_argument('--reference', required=True,
                        help='Reference tRNA FASTA file')
    parser.add_argument('--sample-df', default=None,
                        help='Path to sample_df.xlsx (optional; auto-detects samples if omitted)')
    parser.add_argument('--output-dir', '-o', required=True,
                        help='Output directory for results')
    parser.add_argument('--organism', default='Escherichia coli',
                        help='Organism name for MODOMICS lookup (default: "Escherichia coli")')
    parser.add_argument('--min-coverage', type=int, default=50,
                        help='Minimum read coverage per position (default: 50)')
    parser.add_argument('--n-jobs', type=int, default=4,
                        help='Number of parallel jobs (default: 4)')
    parser.add_argument('--no-modomics', action='store_true',
                        help='Skip MODOMICS API, use fallback CSV only')
    parser.add_argument('--discover-novel', action='store_true',
                        help='Enable novel modification discovery (default: off)')
    parser.add_argument('--csv', action='store_true',
                        help='Also write CSV copies (default: parquet only)')
    parser.set_defaults(func=run_modifications)


def run_modifications(args):
    """Execute the modifications command."""
    from pathlib import Path
    import pandas as pd

    # Import our modules
    from trnaseq.modifications.positional import PositionalExtractor
    from trnaseq.modifications.modomics import MODOMICSAnnotator
    from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer
    from trnaseq.modifications.modification_caller import ModificationCaller

    json_dir = Path(args.json_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"JSON dir:    {json_dir}")
    print(f"Reference:   {args.reference}")
    print(f"Output:      {output_dir}")
    print(f"Organism:    {args.organism}")
    print(f"Min coverage: {args.min_coverage}")
    print(f"Jobs:        {args.n_jobs}")
    print(f"MODOMICS:    {'fallback only' if args.no_modomics else 'API + fallback'}")
    print(f"Novel:       {'enabled' if args.discover_novel else 'disabled'}")
    print()

    # Step 1: Discover sample names
    if args.sample_df:
        sample_df = pd.read_excel(args.sample_df)
        if 'sample_name_unique' in sample_df.columns:
            sample_names = sample_df['sample_name_unique'].tolist()
        else:
            sample_names = sample_df.iloc[:, 0].tolist()
    else:
        # Auto-detect from JSON files
        # .name is e.g. '0p_SWalign.json.bz2' — strip both suffixes + _SWalign
        sample_names = [
            p.name.replace('_SWalign.json.bz2', '')
            for p in json_dir.glob('*_SWalign.json.bz2')
            if not p.name.startswith('common-seqs')
        ]

    print(f"Found {len(sample_names)} samples")

    # Step 2: Extract PSCMs
    print("\nExtracting per-position count matrices...")
    extractor = PositionalExtractor(args.reference)
    all_pscm = extractor.run_parallel(json_dir, sample_names, n_jobs=args.n_jobs)

    # Step 3: Load MODOMICS
    print("\nLoading modification database...")
    annotator = MODOMICSAnnotator(args.organism)
    mods_df = annotator.get_modifications(use_api=not args.no_modomics)
    print(f"  {len(mods_df)} known modification entries loaded")

    # Save MODOMICS annotation
    _save_df(mods_df, output_dir / 'modomics_annotation', args.csv)

    # Step 4: Analyze each sample
    analyzer = RTSignatureAnalyzer(
        min_coverage=args.min_coverage,
        verbose=False
    )
    analyzer.load_reference(args.reference)
    caller = ModificationCaller(organism=args.organism)

    for sample_name, pscm_dict in all_pscm.items():
        print(f"\nProcessing {sample_name} ({len(pscm_dict)} tRNAs)...")

        # Convert to analyzer format
        pscm_dfs = analyzer.load_pscm_from_positional(pscm_dict)

        # Save per-sample profiles
        rt_profile = extractor.compute_rt_profile(pscm_dict)
        mismatch_profile = extractor.compute_mismatch_profile(pscm_dict)
        ac_coverage = extractor.compute_anticodon_coverage(pscm_dict)

        sample_dir = output_dir / sample_name
        sample_dir.mkdir(exist_ok=True)

        _save_df(rt_profile, sample_dir / 'rt_profile', args.csv)
        _save_df(mismatch_profile, sample_dir / 'mismatch_profile', args.csv)
        if not ac_coverage.empty:
            _save_df(ac_coverage, sample_dir / 'anticodon_coverage', args.csv)

        # Save stacked PSCM
        pscm_rows = []
        for trna_name, mat in pscm_dict.items():
            ref_seq = extractor.ref_dict[trna_name]['seq']
            for pos_idx in range(mat.shape[0]):
                pscm_rows.append({
                    'tRNA_name': trna_name,
                    'position': pos_idx + 1,
                    'ref_nt': ref_seq[pos_idx] if pos_idx < len(ref_seq) else 'N',
                    'A': int(mat[pos_idx, 0]),
                    'C': int(mat[pos_idx, 1]),
                    'G': int(mat[pos_idx, 2]),
                    'T': int(mat[pos_idx, 3]),
                    'N': int(mat[pos_idx, 4]),
                    'gap': int(mat[pos_idx, 5]),
                    'coverage': int(mat[pos_idx, 6]),
                    'rt_stop': int(mat[pos_idx, 7]),
                })
        pscm_df_stacked = pd.DataFrame(pscm_rows)
        _save_df(pscm_df_stacked, sample_dir / 'pscm', args.csv)

        # Modification calling
        all_mod_calls = []
        for trna_name, pscm_df in pscm_dfs.items():
            raw_mat = pscm_dict.get(trna_name)
            rt_stops = raw_mat[:, 7] if raw_mat is not None else None

            analysis = analyzer.analyze_trna_with_actual_stops(
                trna_name, pscm_df, rt_stop_counts=rt_stops
            )

            # Annotate with MODOMICS
            annotated = annotator.annotate_signatures(
                analysis['signatures'], trna_name
            )

            ref_seq = analyzer.reference_sequences.get(trna_name, {}).get('seq')
            mod_calls = caller.call_all(
                trna_name, annotated, pscm_df, ref_seq,
                discover_novel=args.discover_novel,
                min_coverage=args.min_coverage,
            )

            if not mod_calls.empty:
                mod_calls['sample'] = sample_name
                all_mod_calls.append(mod_calls)

        if all_mod_calls:
            calls_df = pd.concat(all_mod_calls, ignore_index=True)
            _save_df(calls_df, sample_dir / 'modification_calls', args.csv)
            print(f"  {len(calls_df)} modification calls")
        else:
            print(f"  No modification calls")

    print(f"\nDone! Results saved to {output_dir}/")


def _save_df(df, path_stem, write_csv=False):
    """Save DataFrame as parquet (and optionally CSV)."""
    try:
        df.to_parquet(f'{path_stem}.parquet', index=False)
    except Exception:
        # Fallback to CSV if parquet not available
        df.to_csv(f'{path_stem}.csv', index=False)
        return

    if write_csv:
        df.to_csv(f'{path_stem}.csv', index=False)
