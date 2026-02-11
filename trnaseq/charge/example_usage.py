"""
Example Usage of ChargeQuantifier

This script demonstrates how to use the ChargeQuantifier class
for analyzing tRNA charge levels.

Author: Charge-Specialist
"""

from trnaseq.charge import ChargeQuantifier
import pandas as pd


def basic_usage():
    """
    Basic usage example: Load data and quantify charge.
    """
    print("=" * 70)
    print("Example 1: Basic Usage")
    print("=" * 70)

    # Initialize the quantifier
    quantifier = ChargeQuantifier(
        stats_csv='ALL_stats_aggregate.csv',
        charge_count='count',       # Use read counts for charge calculation
        RPM_count='UMIcount'        # Use UMI counts for RPM calculation
    )

    # Get charge data at transcript level
    charge_df = quantifier.quantify_all(level='transcript')

    print(f"\nLoaded {len(charge_df)} transcript entries")
    print("\nFirst few rows:")
    print(charge_df.head())

    # Export to CSV
    quantifier.export_to_csv(
        output_file='charge_data.csv',
        level='transcript'
    )
    print("\n✓ Exported charge data to charge_data.csv")


def amino_acid_level():
    """
    Example: Analyze charge at amino acid level.
    """
    print("\n" + "=" * 70)
    print("Example 2: Amino Acid Level Analysis")
    print("=" * 70)

    quantifier = ChargeQuantifier(stats_csv='ALL_stats_aggregate.csv')

    # Get charge data at amino acid level
    charge_df = quantifier.quantify_all(level='aa')

    print(f"\nFound {len(charge_df)} amino acid entries")

    # Get summary statistics
    summary = quantifier.get_summary_statistics(level='aa')

    print("\nCharge Statistics by Amino Acid:")
    print(summary[['amino_acid', 'charge_canonical_mean', 'charge_canonical_std']])


def single_trna_analysis():
    """
    Example: Analyze a specific tRNA across samples.
    """
    print("\n" + "=" * 70)
    print("Example 3: Single tRNA Analysis")
    print("=" * 70)

    quantifier = ChargeQuantifier(stats_csv='ALL_stats_aggregate.csv')

    # Analyze Alanine tRNA across all samples
    ala_charge = quantifier.quantify_single_trna(
        trna_id='Ala',
        level='aa'
    )

    print(f"\nAlanine tRNA charge across {len(ala_charge)} samples:")
    print(ala_charge[['sample_name', 'charge_canonical', 'RPM']])


def filtered_analysis():
    """
    Example: Use filtering options.
    """
    print("\n" + "=" * 70)
    print("Example 4: Filtered Analysis")
    print("=" * 70)

    # Exclude alignments with gaps and low scores
    quantifier = ChargeQuantifier(
        stats_csv='ALL_stats_aggregate.csv',
        excl_align_gap=True,        # Exclude alignments with gaps
        excl_09_fmax=True           # Exclude low-quality alignments
    )

    # Get charge data excluding mitochondrial tRNAs
    charge_df = quantifier.quantify_all(
        level='codon',
        include_synthetic=False,    # Exclude synthetic controls
        include_mito=False          # Exclude mitochondrial tRNAs
    )

    print(f"\nFiltered data contains {len(charge_df)} codon entries")
    print("\nCodon-level charge data:")
    print(charge_df[['AA_codon', 'charge_canonical', 'RPM']].head(10))


def compare_samples():
    """
    Example: Compare charge between samples.
    """
    print("\n" + "=" * 70)
    print("Example 5: Compare Samples")
    print("=" * 70)

    quantifier = ChargeQuantifier(stats_csv='ALL_stats_aggregate.csv')

    # Get charge data
    charge_df = quantifier.quantify_all(level='aa')

    # Pivot to compare samples
    pivot_df = charge_df.pivot_table(
        index='amino_acid',
        columns='sample_name',
        values='charge_canonical',
        aggfunc='mean'
    )

    print("\nCharge comparison across samples (first 5 amino acids):")
    print(pivot_df.head())


def main():
    """Run all examples."""
    print("\n" + "=" * 70)
    print("ChargeQuantifier Examples")
    print("=" * 70)

    examples = [
        ("Basic Usage", basic_usage),
        ("Amino Acid Level", amino_acid_level),
        ("Single tRNA Analysis", single_trna_analysis),
        ("Filtered Analysis", filtered_analysis),
        ("Compare Samples", compare_samples)
    ]

    for name, func in examples:
        try:
            func()
        except FileNotFoundError:
            print(f"\n⚠ Skipping '{name}' - CSV file not found")
        except Exception as e:
            print(f"\n⚠ Error in '{name}': {e}")

    print("\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70 + "\n")


if __name__ == '__main__':
    main()
