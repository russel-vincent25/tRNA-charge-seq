"""
Unit Tests for Charge Report + tRNA Source Classification

Run tests with:
    python -m pytest tests/test_charge_report.py -v
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from trnaseq.charge.quantifier import classify_trna_source


class TestClassifyTRNASource:
    """Test classify_trna_source() function."""

    def test_host_default(self):
        """Unmatched tRNAs should return 'host'."""
        prefixes = {'Synthetic_': 'synthetic'}
        assert classify_trna_source('Homo_sapiens_tRNA-Ala-AGC-1-1', prefixes) == 'host'

    def test_synthetic(self):
        """Synthetic prefix should be detected."""
        prefixes = {'Synthetic_': 'synthetic'}
        assert classify_trna_source('Synthetic_tRNA-Ala-AGC-1-1', prefixes) == 'synthetic'

    def test_mutant(self):
        """Mutant prefix should be detected."""
        prefixes = {'Synthetic_': 'synthetic', 'Mutant_itRNA_': 'mutant'}
        assert classify_trna_source('Mutant_itRNA_tRNA-Ala-AGC-1-1', prefixes) == 'mutant'

    def test_first_match_wins(self):
        """First matching prefix should win."""
        prefixes = {'Syn': 'synthetic', 'Synthetic_': 'also_synthetic'}
        assert classify_trna_source('Synthetic_tRNA-Ala-AGC-1-1', prefixes) == 'synthetic'

    def test_empty_prefixes(self):
        """Empty prefix dict should return 'host' for everything."""
        assert classify_trna_source('anything', {}) == 'host'

    def test_custom_category(self):
        """Custom category names should work."""
        prefixes = {'eColi_': 'spike_in'}
        assert classify_trna_source('eColi_tRNA-Lys', prefixes) == 'spike_in'


class TestChargeQuantifierSourceIntegration:
    """Test that tRNA_source column is added to ChargeQuantifier output."""

    @pytest.fixture
    def stats_csv(self, tmp_path):
        """Create a minimal stats CSV with synthetic and host tRNAs."""
        data = {
            'sample_name_unique': ['s1'] * 4,
            'sample_name': ['s1'] * 4,
            'replicate': [1] * 4,
            'barcode': ['BC01'] * 4,
            'species': ['human'] * 4,
            'tRNA_annotation': [
                'Homo_sapiens_tRNA-Ala-AGC-1-1',
                'Homo_sapiens_tRNA-Ala-AGC-1-1',
                'Synthetic_tRNA-Lys-CUU-1-1',
                'Synthetic_tRNA-Lys-CUU-1-1',
            ],
            'tRNA_annotation_len': [73, 73, 70, 70],
            'unique_annotation': [True] * 4,
            '5p_cover': [True] * 4,
            'align_3p_nts': ['CA', 'CC', 'CA', 'CC'],
            'codon': ['GCT', 'GCT', 'AAG', 'AAG'],
            'anticodon': ['AGC', 'AGC', 'CUU', 'CUU'],
            'amino_acid': ['Ala', 'Ala', 'Lys', 'Lys'],
            'count': [100, 50, 30, 10],
        }
        csv_path = tmp_path / 'test_stats.csv'
        pd.DataFrame(data).to_csv(csv_path, index=False)
        return str(csv_path)

    def test_trna_source_column_exists(self, stats_csv):
        from trnaseq.charge.quantifier import ChargeQuantifier
        q = ChargeQuantifier(stats_csv)
        assert 'tRNA_source' in q.stats_df.columns

    def test_trna_source_values(self, stats_csv):
        from trnaseq.charge.quantifier import ChargeQuantifier
        q = ChargeQuantifier(stats_csv)
        sources = q.stats_df['tRNA_source'].unique()
        assert 'host' in sources
        assert 'synthetic' in sources

    def test_syn_ctr_backward_compat(self, stats_csv):
        from trnaseq.charge.quantifier import ChargeQuantifier
        q = ChargeQuantifier(stats_csv)
        assert 'Syn_ctr' in q.stats_df.columns
        syn_rows = q.stats_df[q.stats_df['tRNA_source'] == 'synthetic']
        assert syn_rows['Syn_ctr'].all()

    def test_custom_source_prefixes(self, stats_csv):
        from trnaseq.charge.quantifier import ChargeQuantifier
        q = ChargeQuantifier(
            stats_csv,
            source_prefixes={'Synthetic_': 'spike_in', 'Homo_': 'human'}
        )
        sources = set(q.stats_df['tRNA_source'].unique())
        assert 'spike_in' in sources
        assert 'human' in sources

    def test_trna_source_in_charge_output(self, stats_csv):
        from trnaseq.charge.quantifier import ChargeQuantifier
        q = ChargeQuantifier(stats_csv)
        # transcript-level output should have tRNA_source
        df_tr = q.quantify_all(level='transcript', include_synthetic=True)
        assert 'tRNA_source' in df_tr.columns

    def test_trna_source_in_aa_output(self, stats_csv):
        from trnaseq.charge.quantifier import ChargeQuantifier
        q = ChargeQuantifier(stats_csv)
        df_aa = q.quantify_all(level='aa', include_synthetic=True)
        assert 'tRNA_source' in df_aa.columns


class TestChargeReportGenerator:
    """Test ChargeReportGenerator class."""

    @pytest.fixture
    def charge_data(self):
        """Create minimal charge DataFrames for report generation."""
        np.random.seed(42)
        samples = ['s1', 's1', 's2', 's2']
        sample_names = ['cond1', 'cond1', 'cond2', 'cond2']
        aas = ['Ala', 'Gly', 'Ala', 'Gly']

        charge_tr = pd.DataFrame({
            'sample_name_unique': samples * 2,
            'sample_name': sample_names * 2,
            'tRNA_annotation': ['tRNA-Ala-AGC-1', 'tRNA-Gly-GCC-1'] * 4,
            'tRNA_anno_short': ['Ala-AGC-1', 'Gly-GCC-1'] * 4,
            'amino_acid': aas * 2,
            'charge_canonical': np.random.uniform(30, 90, 8),
            'RPM': np.random.uniform(100, 10000, 8),
            'tRNA_source': ['host'] * 8,
            'Syn_ctr': [False] * 8,
        })

        charge_aa = pd.DataFrame({
            'sample_name_unique': samples,
            'amino_acid': aas,
            'charge_canonical': np.random.uniform(40, 80, 4),
            'tRNA_source': ['host'] * 4,
            'Syn_ctr': [False] * 4,
        })

        charge_summary = pd.DataFrame({
            'amino_acid': ['Ala', 'Gly'],
            'charge_canonical_mean': [65.0, 55.0],
            'level': ['aa', 'aa'],
        })

        sample_df = pd.DataFrame({
            'sample_name_unique': ['s1', 's2'],
            'sample_name': ['cond1', 'cond2'],
        })

        return charge_tr, charge_aa, charge_summary, sample_df

    def test_report_generation(self, charge_data, tmp_path):
        from trnaseq.qc.charge_report import ChargeReportGenerator
        charge_tr, charge_aa, charge_summary, sample_df = charge_data

        gen = ChargeReportGenerator(
            charge_df_transcript=charge_tr,
            charge_df_aa=charge_aa,
            charge_summary=charge_summary,
            sample_df=sample_df,
        )
        out = tmp_path / 'charge_report.html'
        result = gen.generate_html_report(out)
        assert Path(result).exists()
        content = Path(result).read_text()
        assert 'tRNA Charge Dashboard' in content
        assert 'Mean Charge per Sample' in content

    def test_report_panels_present(self, charge_data, tmp_path):
        from trnaseq.qc.charge_report import ChargeReportGenerator
        charge_tr, charge_aa, charge_summary, sample_df = charge_data

        gen = ChargeReportGenerator(
            charge_df_transcript=charge_tr,
            charge_df_aa=charge_aa,
            charge_summary=charge_summary,
            sample_df=sample_df,
        )
        out = tmp_path / 'charge_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Charge by Amino Acid' in content
        assert 'Charge vs Abundance' in content

    def test_empty_data(self, tmp_path):
        from trnaseq.qc.charge_report import ChargeReportGenerator
        gen = ChargeReportGenerator(
            charge_df_transcript=pd.DataFrame(),
            charge_df_aa=pd.DataFrame(),
            charge_summary=pd.DataFrame(),
            sample_df=pd.DataFrame(columns=['sample_name_unique', 'sample_name']),
        )
        out = tmp_path / 'charge_report.html'
        result = gen.generate_html_report(out)
        assert Path(result).exists()

    def test_synthetic_control_panel(self, charge_data, tmp_path):
        """Synthetic control panel renders when prefixes match."""
        from trnaseq.qc.charge_report import ChargeReportGenerator
        charge_tr, charge_aa, charge_summary, sample_df = charge_data

        # Add synthetic rows
        syn_rows = pd.DataFrame({
            'sample_name_unique': ['s1', 's2'],
            'sample_name': ['cond1', 'cond2'],
            'tRNA_annotation': ['Synthetic_tRNA-Lys-CUU-1'] * 2,
            'tRNA_anno_short': ['Lys-CUU-1'] * 2,
            'amino_acid': ['Lys'] * 2,
            'charge_canonical': [2.5, 1.8],
            'RPM': [500, 600],
            'tRNA_source': ['synthetic'] * 2,
            'Syn_ctr': [True] * 2,
        })
        charge_tr_with_syn = pd.concat([charge_tr, syn_rows], ignore_index=True)

        gen = ChargeReportGenerator(
            charge_df_transcript=charge_tr_with_syn,
            charge_df_aa=charge_aa,
            charge_summary=charge_summary,
            sample_df=sample_df,
            source_prefixes={'Synthetic_': 'synthetic'},
        )
        out = tmp_path / 'charge_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Synthetic Control Charge' in content

    def test_synthetic_panel_skipped_without_prefixes(self, charge_data, tmp_path):
        """Synthetic panel not shown when source_prefixes is None."""
        from trnaseq.qc.charge_report import ChargeReportGenerator
        charge_tr, charge_aa, charge_summary, sample_df = charge_data

        gen = ChargeReportGenerator(
            charge_df_transcript=charge_tr,
            charge_df_aa=charge_aa,
            charge_summary=charge_summary,
            sample_df=sample_df,
        )
        out = tmp_path / 'charge_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Synthetic Control' not in content
