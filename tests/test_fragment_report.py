"""
Unit Tests for Fragment Report Generator

Run tests with:
    python -m pytest tests/test_fragment_report.py -v
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path


class TestFragmentReportGenerator:
    """Test FragmentReportGenerator class."""

    @pytest.fixture
    def fragment_data(self):
        """Create minimal fragment DataFrames for report generation."""
        np.random.seed(42)
        samples = ['s1', 's2']

        fragment_counts = pd.DataFrame({
            'sample_name_unique': ['s1', 's1', 's2', 's2'],
            'tRNA_annotation': ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1'] * 2,
            'amino_acid': ['Ala', 'Gly'] * 2,
            'total_reads': [1000, 800, 1200, 600],
            'integrity_score': [0.85, 0.72, 0.90, 0.65],
            # frac_* columns drive the fragment-composition panel; each row sums to 1.
            'frac_full_length': [0.65, 0.55, 0.70, 0.50],
            'frac_rt_dropoff': [0.20, 0.25, 0.18, 0.30],
            'frac_5p_tRF': [0.10, 0.12, 0.08, 0.13],
            'frac_degraded': [0.05, 0.08, 0.04, 0.07],
        })

        rt_dropoff = pd.DataFrame({
            'sample_name_unique': ['s1'] * 10 + ['s2'] * 10,
            'tRNA_annotation': ['tRNA-Ala-AGC-1'] * 10 + ['tRNA-Gly-GCC-1'] * 10,
            'position': list(range(10)) * 2,
            'rt_stop_fraction': np.random.uniform(0, 0.2, 20),
        })

        fragment_lengths = pd.DataFrame({
            'fragment_length': list(range(30, 80)) * 2,
            'read_count': np.random.randint(10, 500, 100),
            'fragment_type': ['full_length'] * 50 + ['rt_dropoff'] * 50,
        })

        fragment_summary = pd.DataFrame({
            'sample_name_unique': samples,
            'N_full_length': [5000, 6000],
            'N_rt_dropoff': [2000, 1500],
            'N_5p_tRF': [500, 300],
            'N_degraded': [200, 150],
            'pct_full_length': [64.9, 75.5],
            'pct_rt_dropoff': [26.0, 18.9],
        })

        sample_df = pd.DataFrame({
            'sample_name_unique': samples,
            'sample_name': ['cond1', 'cond2'],
        })

        return fragment_counts, rt_dropoff, fragment_lengths, fragment_summary, sample_df

    @pytest.fixture
    def coverage_data(self):
        """Create coverage DataFrame for Behrens/needle panel tests."""
        rows = []
        for sample in ['s1', 's2']:
            for aa in ['Ala', 'Gly', 'Lys']:
                max_len = 75
                for pos in range(max_len):
                    rows.append({
                        'sample_name_unique': sample,
                        'amino_acid': aa,
                        'position': pos,
                        'count': max(0, 100 - pos * 2 + np.random.randint(-5, 5)),
                        'max_tRNA_len': max_len,
                    })
        return pd.DataFrame(rows)

    def test_report_generation(self, fragment_data, tmp_path):
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        counts, rt, lengths, summary, sample_df = fragment_data

        gen = FragmentReportGenerator(
            fragment_counts_df=counts,
            rt_dropoff_df=rt,
            fragment_lengths_df=lengths,
            fragment_summary_df=summary,
            sample_df=sample_df,
        )
        out = tmp_path / 'fragment_report.html'
        result = gen.generate_html_report(out)
        assert Path(result).exists()
        content = Path(result).read_text()
        assert 'tRNA Fragment Dashboard' in content

    def test_all_panels_present(self, fragment_data, coverage_data, tmp_path):
        """Every panel in the current dashboard renders when fully fed.

        The RT drop-off profile and fragment length distribution panels were
        deliberately removed in the dashboard overhaul (commit 5c9e122), so
        they are not asserted here.
        """
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        counts, rt, lengths, summary, sample_df = fragment_data

        syn_rows = pd.DataFrame({
            'sample_name_unique': ['s1', 's2'],
            'tRNA_annotation': ['Synthetic_tRNA-Lys-CUU-1'] * 2,
            'amino_acid': ['Lys'] * 2,
            'total_reads': [500, 400],
            'integrity_score': [0.95, 0.93],
            'frac_full_length': [0.90, 0.88],
            'frac_rt_dropoff': [0.06, 0.07],
            'frac_5p_tRF': [0.03, 0.03],
            'frac_degraded': [0.01, 0.02],
        })
        counts_with_syn = pd.concat([counts, syn_rows], ignore_index=True)

        gen = FragmentReportGenerator(
            fragment_counts_df=counts_with_syn,
            rt_dropoff_df=rt,
            fragment_lengths_df=lengths,
            fragment_summary_df=summary,
            sample_df=sample_df,
            coverage_df=coverage_data,
            source_prefixes={'Synthetic_': 'synthetic'},
        )
        out = tmp_path / 'fragment_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Fragment Types' in content
        assert 'Fragment Composition by tRNA' in content
        assert 'Integrity by Amino Acid' in content
        assert 'Behrens Coverage Plot' in content
        assert 'Needle Coverage Plot' in content
        assert 'Synthetic Control Integrity' in content

    def test_empty_data(self, tmp_path):
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        gen = FragmentReportGenerator(
            fragment_counts_df=pd.DataFrame(),
            rt_dropoff_df=pd.DataFrame(),
            fragment_lengths_df=pd.DataFrame(),
            fragment_summary_df=pd.DataFrame(),
            sample_df=pd.DataFrame(columns=['sample_name_unique', 'sample_name']),
        )
        out = tmp_path / 'fragment_report.html'
        result = gen.generate_html_report(out)
        assert Path(result).exists()

    def test_isoacceptor_sort_key(self):
        from trnaseq.qc.fragment_report import _isoacceptor_sort_key
        assert _isoacceptor_sort_key('Prefix-Ala-AGC-1') == ('Ala', 'AGC')
        assert _isoacceptor_sort_key('Prefix-Gly-GCC-2') == ('Gly', 'GCC')
        assert _isoacceptor_sort_key('unknown')[0] == 'unknown'

    def test_behrens_coverage_panel(self, fragment_data, coverage_data, tmp_path):
        """Behrens coverage panel renders when coverage_df is provided."""
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        counts, rt, lengths, summary, sample_df = fragment_data

        gen = FragmentReportGenerator(
            fragment_counts_df=counts,
            rt_dropoff_df=rt,
            fragment_lengths_df=lengths,
            fragment_summary_df=summary,
            sample_df=sample_df,
            coverage_df=coverage_data,
        )
        out = tmp_path / 'fragment_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Behrens Coverage Plot' in content

    def test_needle_coverage_panel(self, fragment_data, coverage_data, tmp_path):
        """Needle coverage panel renders when coverage_df is provided."""
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        counts, rt, lengths, summary, sample_df = fragment_data

        gen = FragmentReportGenerator(
            fragment_counts_df=counts,
            rt_dropoff_df=rt,
            fragment_lengths_df=lengths,
            fragment_summary_df=summary,
            sample_df=sample_df,
            coverage_df=coverage_data,
        )
        out = tmp_path / 'fragment_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Needle Coverage Plot' in content

    def test_coverage_panels_skipped_without_data(self, fragment_data, tmp_path):
        """Coverage panels gracefully skipped when no coverage_df."""
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        counts, rt, lengths, summary, sample_df = fragment_data

        gen = FragmentReportGenerator(
            fragment_counts_df=counts,
            rt_dropoff_df=rt,
            fragment_lengths_df=lengths,
            fragment_summary_df=summary,
            sample_df=sample_df,
        )
        out = tmp_path / 'fragment_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Behrens' not in content
        assert 'Needle' not in content

    def test_synthetic_integrity_panel(self, fragment_data, tmp_path):
        """Synthetic integrity panel renders when source_prefixes match."""
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        counts, rt, lengths, summary, sample_df = fragment_data

        # Add synthetic tRNAs to counts
        syn_rows = pd.DataFrame({
            'sample_name_unique': ['s1', 's2'],
            'tRNA_annotation': ['Synthetic_tRNA-Lys-CUU-1'] * 2,
            'total_reads': [500, 400],
            'integrity_score': [0.95, 0.93],
        })
        counts_with_syn = pd.concat([counts, syn_rows], ignore_index=True)

        gen = FragmentReportGenerator(
            fragment_counts_df=counts_with_syn,
            rt_dropoff_df=rt,
            fragment_lengths_df=lengths,
            fragment_summary_df=summary,
            sample_df=sample_df,
            source_prefixes={'Synthetic_': 'synthetic'},
        )
        out = tmp_path / 'fragment_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Synthetic Control Integrity' in content

    def test_synthetic_panel_skipped_without_prefixes(self, fragment_data, tmp_path):
        """Synthetic panel not shown when source_prefixes is None."""
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        counts, rt, lengths, summary, sample_df = fragment_data

        gen = FragmentReportGenerator(
            fragment_counts_df=counts,
            rt_dropoff_df=rt,
            fragment_lengths_df=lengths,
            fragment_summary_df=summary,
            sample_df=sample_df,
        )
        out = tmp_path / 'fragment_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Synthetic Control' not in content
