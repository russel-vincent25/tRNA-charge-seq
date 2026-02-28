"""
Tests for tRNA modification analysis pipeline.

Uses real test data from projects/example_script_test/data/SWalign/
with masked human tRNA reference (tRNA_database_masked/human/human-tRNAs.fa).
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
REPO_ROOT = Path(__file__).parent.parent
SWALIGN_DIR = REPO_ROOT / 'projects' / 'example_script_test' / 'data' / 'SWalign'
REF_FASTA = REPO_ROOT / 'tRNA_database_masked' / 'human' / 'human-tRNAs.fa'
FALLBACK_DIR = REPO_ROOT / 'trnaseq' / 'modifications' / 'data'
TEST_JSON = SWALIGN_DIR / '0p_SWalign.json.bz2'


# Skip all tests if test data doesn't exist
pytestmark = pytest.mark.skipif(
    not TEST_JSON.exists() or not REF_FASTA.exists(),
    reason="Test data not available"
)


class TestPositionalExtractor:
    """Integration tests for PositionalExtractor with real data."""

    def test_extract_sample_produces_pscm(self):
        """Run extract_sample on 0p_SWalign.json.bz2, verify PSCM output."""
        from trnaseq.modifications.positional import PositionalExtractor

        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)

        # Should have some tRNAs
        assert len(pscm) > 0, "PSCM should contain at least one tRNA"

        # Each entry should be ndarray with shape (ref_len, 8)
        for trna_name, mat in pscm.items():
            assert mat.ndim == 2
            assert mat.shape[1] == 8
            assert mat.shape[0] == ext.ref_dict[trna_name]['seq_len']
            # Coverage (col 6) should be non-negative
            assert np.all(mat[:, 6] >= 0)
            # RT stops (col 7) should be non-negative
            assert np.all(mat[:, 7] >= 0)

    def test_coverage_consistency(self):
        """Coverage column should equal sum of nt counts + gaps."""
        from trnaseq.modifications.positional import PositionalExtractor

        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)

        for trna_name, mat in list(pscm.items())[:5]:  # Check first 5
            nt_sum = mat[:, 0] + mat[:, 1] + mat[:, 2] + mat[:, 3] + mat[:, 4] + mat[:, 5]
            np.testing.assert_array_almost_equal(
                mat[:, 6], nt_sum,
                err_msg=f"Coverage mismatch for {trna_name}"
            )

    def test_rt_profile(self):
        """RT profile should have expected columns and valid values."""
        from trnaseq.modifications.positional import PositionalExtractor

        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)

        rt_df = ext.compute_rt_profile(pscm)

        assert not rt_df.empty
        assert 'tRNA_name' in rt_df.columns
        assert 'position' in rt_df.columns
        assert 'rt_stop_count' in rt_df.columns
        assert 'rt_stop_fraction' in rt_df.columns
        assert 'ref_nt' in rt_df.columns

        # All positions should be >= 1
        assert rt_df['position'].min() >= 1
        # RT stop fraction should be between 0 and 1
        assert rt_df['rt_stop_fraction'].min() >= 0
        assert rt_df['rt_stop_fraction'].max() <= 1.0

    def test_mismatch_profile(self):
        """Mismatch profile should have expected columns."""
        from trnaseq.modifications.positional import PositionalExtractor

        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)

        mm_df = ext.compute_mismatch_profile(pscm)

        assert not mm_df.empty
        assert 'mismatch_rate' in mm_df.columns
        assert 'deletion_rate' in mm_df.columns
        assert 'A_to_C' in mm_df.columns  # Substitution columns

        # Mismatch rate should be between 0 and 1
        assert mm_df['mismatch_rate'].min() >= 0
        assert mm_df['mismatch_rate'].max() <= 1.0

    def test_n_masked_positions_zero_mismatch(self):
        """Positions where ref_nt == 'N' (masked) should have mismatch_rate = 0."""
        from trnaseq.modifications.positional import PositionalExtractor

        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)

        mm_df = ext.compute_mismatch_profile(pscm)

        n_positions = mm_df[mm_df['ref_nt'] == 'N']
        if not n_positions.empty:
            assert (n_positions['mismatch_rate'] == 0).all(), \
                "N-masked positions should have mismatch_rate = 0"

    def test_anticodon_coverage(self):
        """Anticodon coverage should produce valid fractions."""
        from trnaseq.modifications.positional import PositionalExtractor

        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)

        ac_df = ext.compute_anticodon_coverage(pscm)

        if not ac_df.empty:
            assert 'fraction_covering' in ac_df.columns
            assert ac_df['fraction_covering'].min() >= 0
            assert ac_df['fraction_covering'].max() <= 1.0


class TestMODOMICSFallback:
    """Test MODOMICS fallback CSV loading."""

    def test_load_ecoli_fallback(self):
        """Load E. coli fallback CSV and verify structure."""
        from trnaseq.modifications.modomics import MODOMICSAnnotator

        ann = MODOMICSAnnotator('Escherichia coli', fallback_dir=FALLBACK_DIR)
        df = ann.load_fallback()

        assert not df.empty
        assert 'modification_short_name' in df.columns
        assert 'rt_signature_type' in df.columns

        # Should contain well-known E. coli modifications
        mods = set(df['modification_short_name'].values)
        assert 'm1A' in mods or 's4U' in mods or 'm7G' in mods

    def test_load_human_fallback(self):
        """Load human fallback CSV."""
        from trnaseq.modifications.modomics import MODOMICSAnnotator

        ann = MODOMICSAnnotator('Homo sapiens', fallback_dir=FALLBACK_DIR)
        df = ann.load_fallback()

        assert not df.empty
        mods = set(df['modification_short_name'].values)
        assert 'm1A' in mods, "Human fallback should contain m1A"
        assert 'm7G' in mods, "Human fallback should contain m7G"

    def test_load_mouse_fallback(self):
        """Load mouse fallback CSV."""
        from trnaseq.modifications.modomics import MODOMICSAnnotator

        ann = MODOMICSAnnotator('Mus musculus', fallback_dir=FALLBACK_DIR)
        df = ann.load_fallback()

        assert not df.empty

    def test_all_fallback_organisms(self):
        """All 3 organisms should load successfully."""
        from trnaseq.modifications.modomics import MODOMICSAnnotator

        for organism in ['Escherichia coli', 'Homo sapiens', 'Mus musculus']:
            ann = MODOMICSAnnotator(organism, fallback_dir=FALLBACK_DIR)
            df = ann.load_fallback()
            assert not df.empty, f"Fallback for {organism} should not be empty"

            # All entries should have valid rt_signature_type
            valid_types = {'rt_stop', 'mismatch', 'combined', 'silent'}
            actual_types = set(df['rt_signature_type'].values)
            assert actual_types.issubset(valid_types), \
                f"Invalid rt_signature_type values for {organism}: {actual_types - valid_types}"

    def test_invalid_organism_raises(self):
        """Unknown organism should raise FileNotFoundError."""
        from trnaseq.modifications.modomics import MODOMICSAnnotator

        ann = MODOMICSAnnotator('Thermus aquaticus', fallback_dir=FALLBACK_DIR)
        with pytest.raises(FileNotFoundError):
            ann.load_fallback()


class TestModificationCaller:
    """Test modification calling with novel discovery."""

    def test_benjamini_hochberg(self):
        """BH FDR correction should work on simple cases."""
        from trnaseq.modifications.modification_caller import benjamini_hochberg_fdr

        # All significant
        pvals = np.array([0.001, 0.002, 0.003])
        result = benjamini_hochberg_fdr(pvals, alpha=0.05)
        assert result.all()

        # None significant
        pvals = np.array([0.5, 0.6, 0.7])
        result = benjamini_hochberg_fdr(pvals, alpha=0.05)
        assert not result.any()

        # Empty
        result = benjamini_hochberg_fdr(np.array([]), alpha=0.05)
        assert len(result) == 0

    def test_known_profiles_exist(self):
        """MODIFICATION_PROFILES should contain expected entries."""
        from trnaseq.modifications.modification_caller import MODIFICATION_PROFILES

        assert 'm1A' in MODIFICATION_PROFILES
        assert 'm7G' in MODIFICATION_PROFILES
        assert 'm3C' in MODIFICATION_PROFILES
        assert 'pseudouridine' in MODIFICATION_PROFILES

    def test_modification_profile_export(self):
        """ModificationProfile should be importable from __init__."""
        from trnaseq.modifications import ModificationProfile
        assert ModificationProfile is not None


class TestRTSignaturesRefactor:
    """Test RT signatures refactoring."""

    def test_load_pscm_from_positional(self):
        """load_pscm_from_positional should convert ndarray to DataFrame."""
        from trnaseq.modifications.positional import PositionalExtractor
        from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer

        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)

        analyzer = RTSignatureAnalyzer(verbose=False)
        analyzer.load_reference(REF_FASTA)

        pscm_dfs = analyzer.load_pscm_from_positional(pscm)

        assert len(pscm_dfs) > 0
        for trna_name, df in list(pscm_dfs.items())[:3]:
            assert isinstance(df, pd.DataFrame)
            assert 'A' in df.columns
            assert 'C' in df.columns
            assert '-' in df.columns

    def test_calculate_rt_stops_with_actual_counts(self):
        """calculate_rt_stops should accept rt_stop_counts parameter."""
        from trnaseq.modifications.positional import PositionalExtractor
        from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer

        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)

        analyzer = RTSignatureAnalyzer(verbose=False)
        analyzer.load_reference(REF_FASTA)
        pscm_dfs = analyzer.load_pscm_from_positional(pscm)

        # Get first tRNA with data
        trna_name = list(pscm_dfs.keys())[0]
        pscm_df = pscm_dfs[trna_name]
        rt_counts = pscm[trna_name][:, 7]

        rt_df = analyzer.calculate_rt_stops(pscm_df, rt_stop_counts=rt_counts)

        assert 'rt_stop_count' in rt_df.columns
        assert 'rt_stop_pct' in rt_df.columns


class TestIntegration:
    """End-to-end integration tests."""

    def test_full_pipeline_single_sample(self):
        """Run full analysis pipeline on a single sample."""
        from trnaseq.modifications.positional import PositionalExtractor
        from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer
        from trnaseq.modifications.modification_caller import ModificationCaller
        from trnaseq.modifications.modomics import MODOMICSAnnotator

        # Step 1: Extract PSCM
        ext = PositionalExtractor(REF_FASTA)
        pscm = ext.extract_sample(TEST_JSON)
        assert len(pscm) > 0

        # Step 2: Convert to analyzer format
        analyzer = RTSignatureAnalyzer(min_coverage=10, verbose=False)
        analyzer.load_reference(REF_FASTA)
        pscm_dfs = analyzer.load_pscm_from_positional(pscm)

        # Step 3: Analyze first tRNA
        trna_name = list(pscm_dfs.keys())[0]
        pscm_df = pscm_dfs[trna_name]
        rt_counts = pscm[trna_name][:, 7]

        analysis = analyzer.analyze_trna_with_actual_stops(
            trna_name, pscm_df, rt_stop_counts=rt_counts
        )

        assert 'mismatch' in analysis
        assert 'rt_stops' in analysis
        assert 'signatures' in analysis

        # Step 4: Load MODOMICS fallback
        ann = MODOMICSAnnotator('Homo sapiens', fallback_dir=FALLBACK_DIR)
        mods_df = ann.get_modifications(use_api=False)
        assert not mods_df.empty

        # Step 5: Call modifications
        caller = ModificationCaller(organism='human', min_confidence=0.3)
        ref_seq = analyzer.reference_sequences[trna_name]['seq']

        calls = caller.call_all(
            trna_name,
            analysis['signatures'],
            pscm_df,
            ref_seq,
            discover_novel=True,
            min_coverage=10,
        )

        # May or may not have calls depending on data, but should not error
        assert isinstance(calls, pd.DataFrame)

        # Verify fold_change column is present when calls are made
        if not calls.empty:
            assert 'fold_change' in calls.columns
            assert 'background_error_rate' in calls.columns


class TestEdgeCases:
    """Unit tests with synthetic data for edge cases."""

    def test_empty_pscm_dict(self):
        """Empty PSCM dict should produce empty DataFrames."""
        from trnaseq.modifications.positional import PositionalExtractor

        ext = PositionalExtractor(REF_FASTA)

        rt_df = ext.compute_rt_profile({})
        assert rt_df.empty

        mm_df = ext.compute_mismatch_profile({})
        assert mm_df.empty

        ac_df = ext.compute_anticodon_coverage({})
        assert ac_df.empty

    def test_benjamini_hochberg_single_value(self):
        """BH should handle single p-value."""
        from trnaseq.modifications.modification_caller import benjamini_hochberg_fdr

        result = benjamini_hochberg_fdr(np.array([0.001]), alpha=0.05)
        assert len(result) == 1
        assert result[0] == True
