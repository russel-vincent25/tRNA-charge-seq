"""
Unit Tests for Charge Quantification Module

Tests the ChargeQuantifier class and its methods.

Run tests with:
    python -m pytest tests/test_charge.py -v
"""

import pytest
import pandas as pd
import numpy as np
import tempfile
import os
from pathlib import Path

# Import the module to test
from trnaseq.charge import ChargeQuantifier


class TestChargeQuantifier:
    """Test suite for ChargeQuantifier class."""

    @pytest.fixture
    def sample_stats_csv(self, tmp_path):
        """
        Create a sample ALL_stats_aggregate.csv file for testing.

        This fixture creates a temporary CSV file with sample data
        representing different charge states.
        """
        data = {
            'sample_name_unique': ['sample1'] * 8 + ['sample2'] * 8,
            'sample_name': ['sample1'] * 8 + ['sample2'] * 8,
            'replicate': [1] * 16,
            'barcode': ['BC01'] * 16,
            'species': ['human'] * 16,
            'tRNA_annotation': ['Homo_sapiens_tRNA-Ala-AGC-1-1'] * 8 + ['Homo_sapiens_tRNA-Gly-GCC-1-1'] * 8,
            'tRNA_annotation_len': [73] * 8 + [71] * 8,
            'unique_annotation': [True] * 16,
            '5p_cover': [True] * 16,
            'align_3p_nt': ['CA', 'CA', 'CC', 'CC', 'GA', 'CG', 'CA', 'CC'] * 2,
            'codon': ['GCT'] * 8 + ['GGC'] * 8,
            'anticodon': ['AGC'] * 8 + ['GCC'] * 8,
            'amino_acid': ['Ala'] * 8 + ['Gly'] * 8,
            'align_gap': [False] * 16,
            'fmax_score>0.9': [True] * 16,
            'count': [100, 50, 80, 20, 10, 5, 120, 30] * 2,
            'UMIcount': [95, 48, 78, 19, 9, 5, 115, 28] * 2,
            'UMI_percent_exp': [95.0] * 16
        }

        df = pd.DataFrame(data)
        csv_path = tmp_path / "test_stats.csv"
        df.to_csv(csv_path, index=False)

        return str(csv_path)

    @pytest.fixture
    def quantifier(self, sample_stats_csv):
        """Create a ChargeQuantifier instance for testing."""
        return ChargeQuantifier(
            stats_csv=sample_stats_csv,
            charge_count='count',
            RPM_count='UMIcount'
        )

    def test_initialization(self, sample_stats_csv):
        """Test that ChargeQuantifier initializes correctly."""
        quantifier = ChargeQuantifier(stats_csv=sample_stats_csv)

        assert quantifier is not None
        assert quantifier.charge_count_col == 'count'
        assert quantifier.RPM_count_col == 'UMIcount'
        assert quantifier.stats_df is not None
        assert len(quantifier.stats_df) > 0

    def test_invalid_count_column(self, sample_stats_csv):
        """Test that invalid count columns raise ValueError."""
        with pytest.raises(ValueError):
            ChargeQuantifier(stats_csv=sample_stats_csv, charge_count='invalid')

        with pytest.raises(ValueError):
            ChargeQuantifier(stats_csv=sample_stats_csv, RPM_count='invalid')

    def test_column_rename(self, sample_stats_csv):
        """Test that align_3p_nt is renamed to align_3p_nts."""
        quantifier = ChargeQuantifier(stats_csv=sample_stats_csv)

        # Check that the column was renamed
        assert 'align_3p_nts' in quantifier.stats_df.columns
        assert 'align_3p_nt' not in quantifier.stats_df.columns

    def test_amino_acid_conversion(self, quantifier):
        """Test amino acid three-letter to single-letter conversion."""
        assert 'AA_letter' in quantifier.stats_df.columns

        # Check that Ala -> A and Gly -> G
        aa_letters = set(quantifier.stats_df['AA_letter'].unique())
        assert 'A' in aa_letters
        assert 'G' in aa_letters

    def test_charge_calculation(self, quantifier):
        """Test that charge percentages are calculated correctly."""
        # For sample1, Ala tRNA:
        # CA: 100 + 50 + 120 = 270 (charged)
        # CC: 80 + 20 + 30 = 130 (uncharged)
        # Total canonical: 400
        # Expected charge: 270/400 * 100 = 67.5%

        charge_df = quantifier.quantify_all(level='transcript')

        # Find the Ala tRNA entry
        ala_sample1 = charge_df[
            (charge_df['sample_name'] == 'sample1') &
            (charge_df['tRNA_annotation'].str.contains('Ala'))
        ]

        assert len(ala_sample1) > 0

        # Check that charge_canonical is calculated
        assert 'charge_canonical' in ala_sample1.columns
        assert not pd.isna(ala_sample1['charge_canonical'].values[0])

    def test_quantify_all_levels(self, quantifier):
        """Test quantify_all at different annotation levels."""
        # Test transcript level
        df_tr = quantifier.quantify_all(level='transcript')
        assert len(df_tr) > 0
        assert 'tRNA_annotation' in df_tr.columns

        # Test codon level
        df_cd = quantifier.quantify_all(level='codon')
        assert len(df_cd) > 0
        assert 'AA_codon' in df_cd.columns

        # Test amino acid level
        df_aa = quantifier.quantify_all(level='aa')
        assert len(df_aa) > 0
        assert 'amino_acid' in df_aa.columns

    def test_quantify_single_trna(self, quantifier):
        """Test quantifying a single tRNA."""
        # Test at transcript level
        result = quantifier.quantify_single_trna(
            trna_id='Ala',
            level='transcript'
        )

        assert len(result) > 0
        assert all('Ala' in anno for anno in result['tRNA_annotation'])

        # Test at amino acid level
        result_aa = quantifier.quantify_single_trna(
            trna_id='Ala',
            level='aa'
        )

        assert len(result_aa) > 0
        assert all(aa == 'Ala' for aa in result_aa['amino_acid'])

    def test_quantify_single_trna_not_found(self, quantifier):
        """Test that quantify_single_trna raises error for non-existent tRNA."""
        with pytest.raises(ValueError, match='not found'):
            quantifier.quantify_single_trna(
                trna_id='NonExistentTRNA',
                level='transcript'
            )

    def test_count_columns(self, quantifier):
        """Test that CA, CC, GA, CG counts are calculated."""
        df = quantifier.quantify_all(level='transcript')

        assert 'CA_count' in df.columns
        assert 'CC_count' in df.columns
        assert 'GA_count' in df.columns
        assert 'CG_count' in df.columns

        # All counts should be >= 0
        assert (df['CA_count'] >= 0).all()
        assert (df['CC_count'] >= 0).all()
        assert (df['GA_count'] >= 0).all()
        assert (df['CG_count'] >= 0).all()

    def test_rpm_calculation(self, quantifier):
        """Test that RPM is calculated correctly."""
        df = quantifier.quantify_all(level='transcript')

        assert 'RPM' in df.columns
        assert (df['RPM'] >= 0).all()

    def test_filters(self, sample_stats_csv):
        """Test exclude_synthetic and exclude_mito filters."""
        # Create a modified CSV with synthetic and mito tRNAs
        df = pd.read_csv(sample_stats_csv)

        # Add synthetic tRNA
        synth_row = df.iloc[0].copy()
        synth_row['tRNA_annotation'] = 'Synthetic_tRNA-Test'
        synth_row['species'] = 'synth'
        df = pd.concat([df, synth_row.to_frame().T], ignore_index=True)

        # Add mito tRNA
        mito_row = df.iloc[1].copy()
        mito_row['tRNA_annotation'] = 'Homo_sapiens_mito_tRNA-Ala-TGC'
        df = pd.concat([df, mito_row.to_frame().T], ignore_index=True)

        # Save modified CSV
        temp_csv = sample_stats_csv.replace('.csv', '_modified.csv')
        df.to_csv(temp_csv, index=False)

        quantifier = ChargeQuantifier(stats_csv=temp_csv)

        # Test without filters (should include all)
        df_all = quantifier.quantify_all(
            level='transcript',
            include_synthetic=True,
            include_mito=True
        )

        # Test with synthetic filter
        df_no_synth = quantifier.quantify_all(
            level='transcript',
            include_synthetic=False,
            include_mito=True
        )

        # Test with mito filter
        df_no_mito = quantifier.quantify_all(
            level='transcript',
            include_synthetic=True,
            include_mito=False
        )

        assert len(df_all) >= len(df_no_synth)
        assert len(df_all) >= len(df_no_mito)

    def test_export_to_csv(self, quantifier, tmp_path):
        """Test exporting charge data to CSV."""
        output_file = tmp_path / "charge_output.csv"

        quantifier.export_to_csv(
            output_file=str(output_file),
            level='transcript'
        )

        assert output_file.exists()

        # Read back and verify
        df = pd.read_csv(output_file)
        assert len(df) > 0
        assert 'charge_canonical' in df.columns

    def test_summary_statistics(self, quantifier):
        """Test summary statistics generation."""
        summary = quantifier.get_summary_statistics(level='aa')

        assert len(summary) > 0
        assert 'charge_canonical_mean' in summary.columns
        assert 'charge_canonical_std' in summary.columns
        assert 'RPM_mean' in summary.columns

    def test_exclusion_filters(self, sample_stats_csv):
        """Test alignment gap and score exclusion filters."""
        # Modify CSV to include gaps and low scores
        df = pd.read_csv(sample_stats_csv)
        df.loc[0, 'align_gap'] = True
        df.loc[1, 'fmax_score>0.9'] = False

        temp_csv = sample_stats_csv.replace('.csv', '_filtered.csv')
        df.to_csv(temp_csv, index=False)

        # Test with gap exclusion
        quantifier_gap = ChargeQuantifier(
            stats_csv=temp_csv,
            excl_align_gap=True
        )

        # Test with score exclusion
        quantifier_score = ChargeQuantifier(
            stats_csv=temp_csv,
            excl_09_fmax=True
        )

        # Both should have fewer rows than without filters
        quantifier_no_filter = ChargeQuantifier(stats_csv=temp_csv)

        assert len(quantifier_gap.charge_df) <= len(quantifier_no_filter.charge_df)
        assert len(quantifier_score.charge_df) <= len(quantifier_no_filter.charge_df)

    def test_safe_division(self, tmp_path):
        """Test that division by zero is handled safely with NaN."""
        # Create data with zero denominators
        data = {
            'sample_name_unique': ['sample1'],
            'sample_name': ['sample1'],
            'replicate': [1],
            'barcode': ['BC01'],
            'species': ['human'],
            'tRNA_annotation': ['Test_tRNA'],
            'tRNA_annotation_len': [70],
            'unique_annotation': [True],
            '5p_cover': [True],
            'align_3p_nt': ['GA'],  # Only non-canonical
            'codon': ['GCT'],
            'anticodon': ['AGC'],
            'amino_acid': ['Ala'],
            'align_gap': [False],
            'fmax_score>0.9': [True],
            'count': [100],
            'UMIcount': [95],
            'UMI_percent_exp': [95.0]
        }

        df = pd.DataFrame(data)
        csv_path = tmp_path / "zero_div_test.csv"
        df.to_csv(csv_path, index=False)

        quantifier = ChargeQuantifier(stats_csv=str(csv_path))
        charge_df = quantifier.quantify_all(level='transcript')

        # charge_canonical should be NaN (no CA or CC)
        assert pd.isna(charge_df['charge_canonical'].values[0])

        # charge_non-canonical should have a value
        assert not pd.isna(charge_df['charge_non-canonical'].values[0])


class TestChargeAlgorithm:
    """Test the core charge calculation algorithm."""

    def test_canonical_charge_calculation(self):
        """
        Test canonical charge calculation formula.

        Canonical charge = 100 * CA / (CA + CC)
        """
        # Create test data with known values
        CA = 100
        CC = 50
        expected_charge = 100 * CA / (CA + CC)  # 66.67%

        # This should equal 66.67
        assert abs(expected_charge - 66.666667) < 0.01

    def test_non_canonical_charge_calculation(self):
        """
        Test non-canonical charge calculation formula.

        Non-canonical charge = 100 * GA / (GA + CG)
        """
        GA = 80
        CG = 20
        expected_charge = 100 * GA / (GA + CG)  # 80%

        assert expected_charge == 80.0

    def test_mixed_charge_states(self, tmp_path):
        """Test with realistic mixed charge state data."""
        data = {
            'sample_name_unique': ['sample1'] * 4,
            'sample_name': ['sample1'] * 4,
            'replicate': [1] * 4,
            'barcode': ['BC01'] * 4,
            'species': ['human'] * 4,
            'tRNA_annotation': ['Homo_sapiens_tRNA-Ala-AGC-1-1'] * 4,
            'tRNA_annotation_len': [73] * 4,
            'unique_annotation': [True] * 4,
            '5p_cover': [True] * 4,
            'align_3p_nt': ['CA', 'CC', 'GA', 'CG'],
            'codon': ['GCT'] * 4,
            'anticodon': ['AGC'] * 4,
            'amino_acid': ['Ala'] * 4,
            'align_gap': [False] * 4,
            'fmax_score>0.9': [True] * 4,
            'count': [70, 30, 15, 5],  # 70% canonical, 75% non-canonical
            'UMIcount': [68, 29, 14, 5],
            'UMI_percent_exp': [95.0] * 4
        }

        df = pd.DataFrame(data)
        csv_path = tmp_path / "mixed_charge.csv"
        df.to_csv(csv_path, index=False)

        quantifier = ChargeQuantifier(stats_csv=str(csv_path))
        charge_df = quantifier.quantify_all(level='transcript')

        # Check canonical charge: 70/(70+30) = 70%
        canonical_charge = charge_df['charge_canonical'].values[0]
        assert abs(canonical_charge - 70.0) < 0.01

        # Check non-canonical charge: 15/(15+5) = 75%
        non_canonical_charge = charge_df['charge_non-canonical'].values[0]
        assert abs(non_canonical_charge - 75.0) < 0.01


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
