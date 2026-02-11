"""
Example: Complete Modification Analysis Workflow

This script demonstrates the full workflow from RT signature detection
to modification calling with confidence scoring.
"""

from rt_signatures import RTSignatureAnalyzer
from modification_caller import ModificationCaller, MODIFICATION_PROFILES
import pandas as pd


def example_complete_workflow():
    """
    Example: Complete workflow from RT signatures to modification calls
    """
    print("="*70)
    print("Complete Modification Analysis Workflow")
    print("="*70)

    # Step 1: RT Signature Analysis
    print("\n[1/4] Initializing RT Signature Analyzer...")
    analyzer = RTSignatureAnalyzer(
        min_coverage=50,
        mismatch_threshold=0.10,
        rt_stop_threshold=20.0,
        verbose=True
    )

    # Load reference
    print("\n[2/4] Loading reference sequences...")
    analyzer.load_reference('path/to/hg38-tRNAs.fa')

    # Process sample
    print("\n[3/4] Processing sample data...")
    pscm_dict = analyzer.process_sample_from_pipeline(
        stats_csv='path/to/sample1_stats.csv.bz2',
        umi_trimmed_fastq='path/to/sample1_UMI-trimmed.fastq.bz2',
        species='human',
        use_umi_count=True
    )

    # Analyze RT signatures
    print("\n[4/4] Analyzing RT signatures...")
    rt_results = analyzer.analyze_all_trnas()

    # Step 2: Modification Calling
    print("\n" + "="*70)
    print("Calling Modifications")
    print("="*70)

    caller = ModificationCaller(
        organism='human',
        min_confidence=0.5,
        use_position_priors=True,
        statistical_test=True,
        alpha=0.01
    )

    # Get reference sequences for pattern matching
    ref_seqs = {name: info['seq'] for name, info in analyzer.reference_sequences.items()}

    # Call modifications
    print("\nCalling modifications for all tRNAs...")
    modifications = caller.call_modifications_for_all_trnas(
        rt_signature_results=rt_results,
        pscm_dict=pscm_dict,
        reference_sequences=ref_seqs
    )

    # Display results
    print(f"\nFound {len(modifications)} modification sites")
    print("\nTop 20 modification calls by confidence:")
    print(modifications.head(20)[['trna_name', 'position', 'modification',
                                   'confidence', 'mismatch_rate', 'rt_stop_pct']])

    # Filter high-confidence calls
    high_conf = caller.filter_by_confidence(modifications, min_confidence=0.7)
    print(f"\n{len(high_conf)} high-confidence calls (>0.7):")
    print(high_conf[['trna_name', 'position', 'modification', 'confidence']])

    # Summarize by modification type
    print("\n" + "="*70)
    print("Modification Summary")
    print("="*70)
    summary = caller.summarize_modifications(modifications)
    print(summary)

    # Export results
    modifications.to_csv('modification_calls.csv', index=False)
    high_conf.to_csv('high_confidence_modifications.csv', index=False)
    summary.to_csv('modification_summary.csv', index=False)

    print("\n✅ Results exported to CSV files")


def example_single_trna_analysis():
    """
    Example: Detailed analysis of a single tRNA
    """
    print("\n" + "="*70)
    print("Single tRNA Analysis - tRNA-Thr with known m1A at position 58")
    print("="*70)

    # Initialize (abbreviated)
    analyzer = RTSignatureAnalyzer()
    analyzer.load_reference('path/to/hg38-tRNAs.fa')
    # ... process sample ...

    # Analyze specific tRNA
    trna_name = 'tRNA-Thr-AGT-1-1'
    pscm_df = analyzer.pscm_data[trna_name]
    results = analyzer.analyze_trna(trna_name, pscm_df)

    # Call modifications
    caller = ModificationCaller(organism='human')
    ref_seq = analyzer.reference_sequences[trna_name]['seq']

    modifications = caller.call_modifications_for_trna(
        trna_name,
        results['signatures'],
        pscm_df,
        ref_seq
    )

    print(f"\nModifications called in {trna_name}:")
    print(modifications)

    # Inspect position 58 (known m1A site)
    pos58_calls = modifications[modifications['position'] == 58]
    if not pos58_calls.empty:
        print("\n✅ Position 58 (expected m1A):")
        print(pos58_calls)

        # Get detailed mutation patterns
        patterns = analyzer.extract_mutation_patterns(pscm_df, 58, ref_seq)
        print("\nMutation patterns at position 58:")
        for pattern, count in sorted(patterns.items(), key=lambda x: x[1], reverse=True):
            total = pscm_df.iloc[57].sum()
            pct = (count / total * 100) if total > 0 else 0
            print(f"  {pattern:15s}: {count:8.1f} ({pct:5.1f}%)")


def example_modification_profiles():
    """
    Example: Explore available modification profiles
    """
    print("\n" + "="*70)
    print("Available Modification Profiles")
    print("="*70)

    for name, profile in MODIFICATION_PROFILES.items():
        print(f"\n{profile.name} ({profile.full_name})")
        print(f"  Typical positions: {profile.typical_positions}")
        print(f"  Signature type: {profile.signature_type}")
        print(f"  Pattern: {profile.mismatch_pattern}")
        print(f"  Minimum rate: {profile.min_rate:.2%}")
        if profile.rt_stop_required:
            print(f"  RT stop required: >{profile.min_rt_stop_pct}%")


def example_compare_samples():
    """
    Example: Compare modifications between samples
    """
    print("\n" + "="*70)
    print("Compare Modifications Between Samples")
    print("="*70)

    # Process multiple samples
    samples = ['control_rep1', 'control_rep2', 'treated_rep1', 'treated_rep2']
    all_modifications = []

    analyzer = RTSignatureAnalyzer()
    analyzer.load_reference('path/to/hg38-tRNAs.fa')
    caller = ModificationCaller(organism='human', min_confidence=0.6)

    for sample_name in samples:
        print(f"\nProcessing {sample_name}...")

        # Process sample
        pscm_dict = analyzer.process_sample_from_pipeline(
            stats_csv=f'path/to/{sample_name}_stats.csv.bz2',
            umi_trimmed_fastq=f'path/to/{sample_name}_UMI-trimmed.fastq.bz2',
            species='human'
        )

        rt_results = analyzer.analyze_all_trnas()

        # Call modifications
        ref_seqs = {n: i['seq'] for n, i in analyzer.reference_sequences.items()}
        mods = caller.call_modifications_for_all_trnas(
            rt_results, pscm_dict, ref_seqs
        )

        mods['sample'] = sample_name
        all_modifications.append(mods)

    # Combine and analyze
    combined = pd.concat(all_modifications, ignore_index=True)

    # Find modifications present in treated but not control
    control_sites = set()
    treated_sites = set()

    for _, row in combined.iterrows():
        site = (row['trna_name'], row['position'], row['modification'])
        if 'control' in row['sample']:
            control_sites.add(site)
        elif 'treated' in row['sample']:
            treated_sites.add(site)

    novel_in_treated = treated_sites - control_sites
    print(f"\nNovel modifications in treated samples: {len(novel_in_treated)}")
    for trna, pos, mod in sorted(novel_in_treated):
        print(f"  {trna} position {pos}: {mod}")

    # Export
    combined.to_csv('multi_sample_modifications.csv', index=False)


def example_validation_against_known_sites():
    """
    Example: Validate modification calls against known modification sites
    """
    print("\n" + "="*70)
    print("Validation Against Known Modification Sites")
    print("="*70)

    # Known modifications from literature or MODOMICS
    known_modifications = {
        'tRNA-Thr-AGT-1-1': [(58, 'm1A'), (46, 'm7G')],
        'tRNA-Ala-AGC-1-1': [(58, 'm1A'), (32, 'Ψ')],
        'tRNA-Lys-TTT-1-1': [(58, 'm1A'), (34, 'i6A')],
    }

    # Run analysis (abbreviated)
    analyzer = RTSignatureAnalyzer()
    analyzer.load_reference('path/to/hg38-tRNAs.fa')
    # ... process sample ...
    rt_results = analyzer.analyze_all_trnas()

    caller = ModificationCaller(organism='human')
    ref_seqs = {n: i['seq'] for n, i in analyzer.reference_sequences.items()}
    calls = caller.call_modifications_for_all_trnas(
        rt_results, analyzer.pscm_data, ref_seqs
    )

    # Validate
    print("\nValidation Results:")
    true_positives = 0
    false_negatives = 0

    for trna_name, expected_mods in known_modifications.items():
        print(f"\n{trna_name}:")
        trna_calls = calls[calls['trna_name'] == trna_name]

        for pos, mod_type in expected_mods:
            found = trna_calls[
                (trna_calls['position'] == pos) &
                (trna_calls['modification'] == mod_type)
            ]

            if not found.empty:
                conf = found.iloc[0]['confidence']
                print(f"  ✅ {mod_type} at {pos}: Found (confidence={conf:.2f})")
                true_positives += 1
            else:
                print(f"  ❌ {mod_type} at {pos}: MISSED")
                false_negatives += 1

    sensitivity = true_positives / (true_positives + false_negatives)
    print(f"\nSensitivity: {sensitivity:.1%} ({true_positives}/{true_positives + false_negatives})")


if __name__ == '__main__':
    print("Modification Calling Examples")
    print("=" * 70)
    print("\nNote: Update file paths before running")
    print()

    # Uncomment to run examples:
    # example_complete_workflow()
    # example_single_trna_analysis()
    # example_modification_profiles()
    # example_compare_samples()
    # example_validation_against_known_sites()

    print("\nTo explore modification profiles:")
    example_modification_profiles()
