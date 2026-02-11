"""
Example usage of the RT Signature Analysis module

This script demonstrates how to use the new modular RT signature analyzer
to process tRNA-seq data and identify modification sites.
"""

from rt_signatures import RTSignatureAnalyzer, analyze_rt_signatures


def example_basic_usage():
    """
    Example 1: Basic RT signature analysis from existing pipeline output
    """
    print("="*60)
    print("Example 1: Basic RT signature analysis")
    print("="*60)

    # Initialize analyzer
    analyzer = RTSignatureAnalyzer(
        min_coverage=50,
        mismatch_threshold=0.10,
        rt_stop_threshold=20.0,
        verbose=True
    )

    # Load reference sequences
    analyzer.load_reference('path/to/hg38-tRNAs.fa')

    # Process a sample from the existing pipeline
    pscm_dict = analyzer.process_sample_from_pipeline(
        stats_csv='path/to/sample1_stats.csv.bz2',
        umi_trimmed_fastq='path/to/sample1_UMI-trimmed.fastq.bz2',
        species='human',
        use_umi_count=True,
        unique_anno=True
    )

    # Analyze all tRNAs
    results = analyzer.analyze_all_trnas()

    # Inspect results for a specific tRNA
    trna_name = 'tRNA-Ala-AGC-1-1'
    if trna_name in results:
        print(f"\nAnalysis results for {trna_name}:")
        print("\nPositions with RT signatures:")
        sig_df = results[trna_name]['signatures']
        print(sig_df[sig_df['has_signature']])

        print("\nTop mismatch positions:")
        mismatch_df = results[trna_name]['mismatch']
        print(mismatch_df.nlargest(10, 'mismatch_rate'))

        print("\nRT stop positions:")
        rt_stop_df = results[trna_name]['rt_stops']
        print(rt_stop_df[rt_stop_df['rt_stop_pct'] > 20])


def example_analyze_single_trna():
    """
    Example 2: Analyze a single tRNA with custom PSCM
    """
    print("\n" + "="*60)
    print("Example 2: Single tRNA analysis with custom PSCM")
    print("="*60)

    import pandas as pd
    import numpy as np

    # Initialize analyzer
    analyzer = RTSignatureAnalyzer(min_coverage=30)
    analyzer.load_reference('path/to/reference.fa')

    # Load pre-computed PSCM (from pickle or CSV)
    # This example shows loading from CSV
    pscm_df = pd.read_csv('path/to/trna_pscm.csv', index_col=0)

    # Analyze single tRNA
    trna_name = 'tRNA-Thr-AGT-1-1'
    results = analyzer.analyze_trna(trna_name, pscm_df)

    # Extract mutation patterns at position 58 (common m1A site)
    patterns = analyzer.extract_mutation_patterns(
        pscm_df,
        position=58,
        reference_seq=analyzer.reference_sequences[trna_name]['seq']
    )

    print(f"\nMutation patterns at position 58:")
    for pattern, count in sorted(patterns.items(), key=lambda x: x[1], reverse=True):
        print(f"  {pattern}: {count:.1f}")


def example_quick_init():
    """
    Example 3: Quick initialization and analysis
    """
    print("\n" + "="*60)
    print("Example 3: Quick initialization")
    print("="*60)

    # Quick initialization with reference loading
    analyzer = analyze_rt_signatures(
        reference_fasta='path/to/hg38-tRNAs.fa',
        min_coverage=50,
        mismatch_threshold=0.10
    )

    print(f"Loaded {len(analyzer.reference_sequences)} tRNA references")


def example_export_results():
    """
    Example 4: Export RT signature results to CSV
    """
    print("\n" + "="*60)
    print("Example 4: Export results to CSV")
    print("="*60)

    # Run analysis (abbreviated)
    analyzer = RTSignatureAnalyzer()
    analyzer.load_reference('path/to/reference.fa')
    # ... process samples ...
    results = analyzer.analyze_all_trnas()

    # Export all signatures to CSV
    import pandas as pd

    all_signatures = []
    for trna_name, trna_results in results.items():
        sig_df = trna_results['signatures'].copy()
        sig_df['trna_name'] = trna_name
        all_signatures.append(sig_df)

    combined_df = pd.concat(all_signatures, ignore_index=True)
    combined_df.to_csv('rt_signatures_all_trnas.csv', index=False)
    print(f"\nExported signatures for {len(results)} tRNAs")

    # Export high-confidence modification sites only
    high_conf = combined_df[
        (combined_df['mismatch_rate'] > 0.15) |
        (combined_df['rt_stop_pct'] > 30)
    ]
    high_conf.to_csv('high_confidence_modifications.csv', index=False)
    print(f"Found {len(high_conf)} high-confidence modification sites")


if __name__ == '__main__':
    print("RT Signature Analysis Examples")
    print("=" * 60)
    print("\nNote: Update file paths before running")
    print()

    # Uncomment to run examples:
    # example_basic_usage()
    # example_analyze_single_trna()
    # example_quick_init()
    # example_export_results()

    print("\nExamples complete!")
