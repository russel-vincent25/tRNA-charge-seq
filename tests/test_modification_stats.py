"""
Tests for enhanced modification statistics:
- Background error-rate estimation
- Replicate aggregation (Fisher's method)
- Fold-change computation
- ModificationCaller with custom background
- Modification QC report generation
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers to build synthetic PSCM data
# ---------------------------------------------------------------------------

def _make_pscm(ref_seq, coverage=500, error_rate=0.005):
    """Build a clean PSCM ndarray for a given reference sequence.

    Columns: A(0) C(1) G(2) T(3) N(4) gap(5) coverage(6) rt_stop(7)
    """
    _NT = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    n = len(ref_seq)
    mat = np.zeros((n, 8), dtype=np.float64)
    for i, nt in enumerate(ref_seq.upper()):
        if nt in _NT:
            correct = int(coverage * (1 - error_rate))
            mismatches = coverage - correct
            mat[i, _NT[nt]] = correct
            # Spread mismatches evenly across other nts
            others = [v for k, v in _NT.items() if k != nt]
            per_other = mismatches // len(others)
            for o in others:
                mat[i, o] = per_other
            mat[i, 6] = coverage
        else:
            mat[i, 6] = coverage
            mat[i, 4] = coverage  # N count
    return mat


def _make_ref_dict(names_seqs):
    return {n: {'seq': s, 'seq_len': len(s)} for n, s in names_seqs}


# ===========================================================================
# TestBackgroundEstimation
# ===========================================================================

class TestBackgroundEstimation:

    def test_synthetic_detection(self):
        """Synthetic tRNAs should be used for background."""
        from trnaseq.modifications.modification_caller import estimate_background_error_rate

        ref = 'ACGTACGT'
        pscm = {
            'Synthetic_spike1': _make_pscm(ref, coverage=1000, error_rate=0.003),
            'tRNA-Ala-AGC-1-1': _make_pscm(ref, coverage=1000, error_rate=0.20),
        }
        ref_dict = _make_ref_dict([
            ('Synthetic_spike1', ref),
            ('tRNA-Ala-AGC-1-1', ref),
        ])
        rate, source = estimate_background_error_rate(pscm, ref_dict)
        assert source == 'synthetic'
        # Should be close to 0.003, not 0.20
        assert rate < 0.01

    def test_fallback_to_q25(self):
        """Without synthetic tRNAs, should fall back to 25th percentile."""
        from trnaseq.modifications.modification_caller import estimate_background_error_rate

        ref = 'ACGTACGT'
        pscm = {
            'tRNA-Ala-AGC-1-1': _make_pscm(ref, coverage=1000, error_rate=0.005),
        }
        ref_dict = _make_ref_dict([('tRNA-Ala-AGC-1-1', ref)])
        rate, source = estimate_background_error_rate(pscm, ref_dict)
        assert source == 'empirical_q25'
        assert rate >= 0.001  # floor

    def test_empty_pscm_returns_floor(self):
        """Empty input should return floor value."""
        from trnaseq.modifications.modification_caller import estimate_background_error_rate

        rate, source = estimate_background_error_rate({}, {})
        assert rate == 0.001
        assert source == 'empirical_q25'

    def test_coverage_weighting(self):
        """High-coverage positions should dominate the estimate."""
        from trnaseq.modifications.modification_caller import estimate_background_error_rate

        ref = 'ACGT'
        # One synthetic with high cov + low error, one with low cov + high error
        mat = np.zeros((4, 8), dtype=np.float64)
        # Position 0: A, 10000 reads, 0.1% error
        mat[0, 0] = 9990  # A
        mat[0, 1] = 5; mat[0, 2] = 3; mat[0, 3] = 2
        mat[0, 6] = 10000
        # Position 1: C, 100 reads, 10% error
        mat[1, 1] = 90; mat[1, 0] = 5; mat[1, 2] = 3; mat[1, 3] = 2
        mat[1, 6] = 100
        # Position 2: G, 10000 reads, 0.2% error
        mat[2, 2] = 9980; mat[2, 0] = 10; mat[2, 1] = 5; mat[2, 3] = 5
        mat[2, 6] = 10000
        # Position 3: T, 100 reads, 10% error
        mat[3, 3] = 90; mat[3, 0] = 5; mat[3, 1] = 3; mat[3, 2] = 2
        mat[3, 6] = 100

        pscm = {'Synthetic_test': mat}
        ref_dict = _make_ref_dict([('Synthetic_test', ref)])
        rate, source = estimate_background_error_rate(pscm, ref_dict)
        assert source == 'synthetic'
        # Weighted by coverage: dominated by high-cov low-error positions
        assert rate < 0.02

    def test_min_coverage_filtering(self):
        """Positions below min_coverage should be ignored."""
        from trnaseq.modifications.modification_caller import estimate_background_error_rate

        ref = 'AC'
        mat = np.zeros((2, 8), dtype=np.float64)
        # Position 0: A, 200 reads, clean
        mat[0, 0] = 198; mat[0, 1] = 1; mat[0, 2] = 1
        mat[0, 6] = 200
        # Position 1: C, 10 reads (below min_coverage=50), very noisy
        mat[1, 1] = 5; mat[1, 0] = 5
        mat[1, 6] = 10

        pscm = {'Synthetic_filt': mat}
        ref_dict = _make_ref_dict([('Synthetic_filt', ref)])
        rate, source = estimate_background_error_rate(
            pscm, ref_dict, min_coverage=50
        )
        assert source == 'synthetic'
        # Only position 0 counted: ~2/200 = 0.01
        assert rate < 0.02


# ===========================================================================
# TestReplicateAggregator
# ===========================================================================

class TestReplicateAggregator:

    def _make_calls(self, trna='tRNA-Ala-AGC-1-1', pos=58, mod='m1A',
                    mm_rate=0.25, fc=25.0, cov=500, conf=0.8, pval=1e-10):
        return pd.DataFrame([{
            'trna_name': trna, 'position': pos, 'modification': mod,
            'mismatch_rate': mm_rate, 'fold_change': fc, 'coverage': cov,
            'confidence': conf, 'pvalue': pval,
        }])

    def test_basic_consensus_3_of_4(self):
        """3 of 4 replicates should reach consensus with min_replicates=3."""
        from trnaseq.modifications.modification_caller import ReplicateAggregator

        per_sample = {
            'S1': self._make_calls(),
            'S2': self._make_calls(),
            'S3': self._make_calls(),
            'S4': pd.DataFrame(),  # 4th replicate has no call
        }
        groups = {'CondA': ['S1', 'S2', 'S3', 'S4']}

        agg = ReplicateAggregator(min_replicates=3, alpha=0.01)
        result = agg.aggregate(per_sample, groups)

        assert len(result) == 1
        assert result.iloc[0]['n_replicates_detected'] == 3
        assert result.iloc[0]['n_replicates_total'] == 4
        assert result.iloc[0]['consensus_call'] == True

    def test_below_threshold_2_of_4(self):
        """2 of 4 replicates should NOT reach consensus with min_replicates=3."""
        from trnaseq.modifications.modification_caller import ReplicateAggregator

        per_sample = {
            'S1': self._make_calls(),
            'S2': self._make_calls(),
            'S3': pd.DataFrame(),
            'S4': pd.DataFrame(),
        }
        groups = {'CondA': ['S1', 'S2', 'S3', 'S4']}

        agg = ReplicateAggregator(min_replicates=3, alpha=0.01)
        result = agg.aggregate(per_sample, groups)

        assert len(result) == 1
        assert result.iloc[0]['n_replicates_detected'] == 2
        assert result.iloc[0]['consensus_call'] == False

    def test_fisher_combined_pvalue(self):
        """Fisher combined p should be < individual p-values."""
        from trnaseq.modifications.modification_caller import ReplicateAggregator

        per_sample = {
            'S1': self._make_calls(pval=0.005),
            'S2': self._make_calls(pval=0.005),
            'S3': self._make_calls(pval=0.005),
        }
        groups = {'CondA': ['S1', 'S2', 'S3']}

        agg = ReplicateAggregator(min_replicates=3, alpha=0.01)
        result = agg.aggregate(per_sample, groups)

        assert result.iloc[0]['fisher_combined_pvalue'] < 0.005
        assert result.iloc[0]['fisher_significant'] == True

    def test_single_replicate(self):
        """Single replicate should still produce a row."""
        from trnaseq.modifications.modification_caller import ReplicateAggregator

        per_sample = {'S1': self._make_calls()}
        groups = {'CondA': ['S1']}

        agg = ReplicateAggregator(min_replicates=1, alpha=0.05)
        result = agg.aggregate(per_sample, groups)

        assert len(result) == 1
        assert result.iloc[0]['n_replicates_detected'] == 1

    def test_empty_calls(self):
        """All-empty calls should produce empty DataFrame."""
        from trnaseq.modifications.modification_caller import ReplicateAggregator

        per_sample = {'S1': pd.DataFrame(), 'S2': pd.DataFrame()}
        groups = {'CondA': ['S1', 'S2']}

        agg = ReplicateAggregator()
        result = agg.aggregate(per_sample, groups)
        assert result.empty

    def test_nan_pvalue_handling(self):
        """NaN p-values should not cause crash."""
        from trnaseq.modifications.modification_caller import ReplicateAggregator

        calls = self._make_calls(pval=np.nan)
        per_sample = {'S1': calls, 'S2': calls, 'S3': calls}
        groups = {'CondA': ['S1', 'S2', 'S3']}

        agg = ReplicateAggregator(min_replicates=3, alpha=0.01)
        result = agg.aggregate(per_sample, groups)

        assert len(result) == 1
        # NaN p-values: Fisher can't be computed, so consensus_call=False
        assert result.iloc[0]['consensus_call'] == False


# ===========================================================================
# TestFoldChange
# ===========================================================================

class TestFoldChange:

    def test_fold_change_computed(self):
        """call_all should produce fold_change column."""
        from trnaseq.modifications.modification_caller import ModificationCaller

        caller = ModificationCaller(
            organism='human', min_confidence=0.0,
            background_error_rate=0.01, statistical_test=False,
        )

        # Build minimal signature DataFrame
        sig = pd.DataFrame([{
            'position': 58,
            'has_signature': True,
            'mismatch_rate': 0.25,
            'rt_stop_pct': 30.0,
            'gap_rate': 0.0,
            'coverage': 1000,
        }])

        calls = caller.call_all('tRNA-Ala-AGC-1-1', sig, discover_novel=False)
        if not calls.empty:
            assert 'fold_change' in calls.columns
            assert 'background_error_rate' in calls.columns
            # 0.25 / 0.01 = 25.0
            assert calls.iloc[0]['fold_change'] == pytest.approx(25.0, rel=0.1)

    def test_zero_background_gives_nan(self):
        """Zero background should produce NaN fold_change."""
        from trnaseq.modifications.modification_caller import ModificationCaller

        caller = ModificationCaller(
            organism='human', min_confidence=0.0,
            background_error_rate=0.0, statistical_test=False,
        )

        sig = pd.DataFrame([{
            'position': 58,
            'has_signature': True,
            'mismatch_rate': 0.25,
            'rt_stop_pct': 30.0,
            'gap_rate': 0.0,
            'coverage': 1000,
        }])

        calls = caller.call_all('tRNA-Ala-AGC-1-1', sig, discover_novel=False)
        if not calls.empty and 'fold_change' in calls.columns:
            assert np.isnan(calls.iloc[0]['fold_change'])


# ===========================================================================
# TestModificationCallerBackground
# ===========================================================================

class TestModificationCallerBackground:

    def test_custom_bg_rate_used(self):
        """Custom background rate should be used in binomial test."""
        from trnaseq.modifications.modification_caller import ModificationCaller

        caller = ModificationCaller(background_error_rate=0.05)
        # With 5% background, 6% mismatch at 100 reads should not be significant
        pval = caller.perform_statistical_test(100, 6)
        assert pval > 0.05  # Not significant against 5% background

    def test_default_backward_compat(self):
        """Default background_error_rate should be 0.01."""
        from trnaseq.modifications.modification_caller import ModificationCaller

        caller = ModificationCaller()
        assert caller.background_error_rate == 0.01

    def test_explicit_override(self):
        """Explicit expected_error_rate should override instance default."""
        from trnaseq.modifications.modification_caller import ModificationCaller

        caller = ModificationCaller(background_error_rate=0.05)
        # Override with explicit 0.001
        pval = caller.perform_statistical_test(1000, 50, expected_error_rate=0.001)
        # 5% mismatch vs 0.1% background should be highly significant
        assert pval < 0.001


# ===========================================================================
# TestModificationReport
# ===========================================================================

class TestModificationReport:

    def test_html_generated(self, tmp_path):
        """Report should generate valid HTML file."""
        from trnaseq.qc.modification_report import ModificationReportGenerator

        calls = pd.DataFrame([{
            'trna_name': 'tRNA-Ala-AGC-1-1', 'position': 58,
            'modification': 'm1A', 'mismatch_rate': 0.25,
            'fold_change': 25.0, 'coverage': 500, 'confidence': 0.9,
            'pvalue': 1e-10,
        }])

        gen = ModificationReportGenerator(
            per_sample_calls={'S1': calls, 'S2': calls},
        )
        out = gen.generate_html_report(tmp_path / 'report.html')
        assert Path(out).exists()
        content = Path(out).read_text()
        assert '<html' in content
        assert 'Modification' in content

    def test_empty_data_produces_valid_html(self, tmp_path):
        """Empty data should still produce a valid HTML file."""
        from trnaseq.qc.modification_report import ModificationReportGenerator

        gen = ModificationReportGenerator(per_sample_calls={})
        out = gen.generate_html_report(tmp_path / 'report.html')
        assert Path(out).exists()
        content = Path(out).read_text()
        assert '<html' in content

    def test_report_with_aggregation(self, tmp_path):
        """Report with aggregated + consensus data should not crash."""
        from trnaseq.qc.modification_report import ModificationReportGenerator

        calls = pd.DataFrame([{
            'trna_name': 'tRNA-Ala-AGC-1-1', 'position': 58,
            'modification': 'm1A', 'mismatch_rate': 0.25,
            'fold_change': 25.0, 'coverage': 500, 'confidence': 0.9,
            'pvalue': 1e-10,
        }])

        agg = pd.DataFrame([{
            'sample_name': 'CondA', 'trna_name': 'tRNA-Ala-AGC-1-1',
            'position': 58, 'modification': 'm1A',
            'n_replicates_detected': 3, 'n_replicates_total': 4,
            'fisher_combined_pvalue': 1e-20, 'fisher_significant': True,
            'mean_mismatch_rate': 0.25, 'mean_fold_change': 25.0,
            'mean_coverage': 500.0, 'mean_confidence': 0.9,
            'consensus_call': True,
        }])

        gen = ModificationReportGenerator(
            per_sample_calls={'S1': calls, 'S2': calls},
            aggregated_calls=agg,
            consensus_calls=agg,
            replicate_groups={'CondA': ['S1', 'S2']},
            ref_dict={'tRNA-Ala-AGC-1-1': {'seq': 'A' * 76, 'seq_len': 76}},
        )
        out = gen.generate_html_report(tmp_path / 'report.html')
        assert Path(out).exists()
