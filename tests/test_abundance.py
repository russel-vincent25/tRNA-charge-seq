"""
Unit Tests for Abundance Module + Report

Run tests with:
    python -m pytest tests/test_abundance.py -v
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path


class TestDifferentialAbundanceCountMatrix:
    """Test DifferentialAbundance count matrix building."""

    @pytest.fixture
    def stats_csv(self, tmp_path):
        """Create a minimal stats CSV."""
        data = {
            'sample_name_unique': ['s1'] * 3 + ['s2'] * 3 + ['s3'] * 3,
            'tRNA_annotation': ['tRNA-Ala', 'tRNA-Gly', 'tRNA-Lys'] * 3,
            'amino_acid': ['Ala', 'Gly', 'Lys'] * 3,
            'codon': ['GCT', 'GGC', 'AAG'] * 3,
            'count': [100, 200, 50, 120, 180, 60, 90, 210, 45],
        }
        csv_path = tmp_path / 'stats.csv'
        pd.DataFrame(data).to_csv(csv_path, index=False)
        return str(csv_path)

    @pytest.fixture
    def sample_df(self):
        return pd.DataFrame({
            'sample_name_unique': ['s1', 's2', 's3'],
            'sample_name': ['WT', 'WT', 'KO'],
        })

    def test_count_matrix_shape(self, stats_csv, sample_df):
        from trnaseq.abundance import DifferentialAbundance
        da = DifferentialAbundance(stats_csv, sample_df, level='aa')
        assert da.count_matrix.shape == (3, 3)  # 3 samples x 3 amino acids

    def test_count_matrix_values(self, stats_csv, sample_df):
        from trnaseq.abundance import DifferentialAbundance
        da = DifferentialAbundance(stats_csv, sample_df, level='aa')
        assert da.count_matrix.loc['s1', 'Ala'] == 100
        assert da.count_matrix.loc['s2', 'Gly'] == 180

    def test_control_group_auto_detect(self, stats_csv, sample_df):
        from trnaseq.abundance import DifferentialAbundance
        da = DifferentialAbundance(stats_csv, sample_df, level='aa')
        assert da.control_group == 'KO'  # first alphabetically

    def test_control_group_explicit(self, stats_csv, sample_df):
        from trnaseq.abundance import DifferentialAbundance
        da = DifferentialAbundance(stats_csv, sample_df, level='aa',
                                   control_group='WT')
        assert da.control_group == 'WT'

    def test_invalid_control_group(self, stats_csv, sample_df):
        from trnaseq.abundance import DifferentialAbundance
        with pytest.raises(ValueError, match="not in groups"):
            DifferentialAbundance(stats_csv, sample_df, level='aa',
                                  control_group='NONEXISTENT')

    def test_invalid_level(self, stats_csv, sample_df):
        from trnaseq.abundance import DifferentialAbundance
        with pytest.raises(ValueError, match="level must be"):
            DifferentialAbundance(stats_csv, sample_df, level='invalid')

    def test_transcript_level(self, stats_csv, sample_df):
        from trnaseq.abundance import DifferentialAbundance
        da = DifferentialAbundance(stats_csv, sample_df, level='transcript')
        assert da.count_matrix.shape[1] == 3  # 3 tRNA annotations


class TestAbundanceReportGenerator:
    """Test AbundanceReportGenerator class."""

    @pytest.fixture
    def report_data(self):
        """Create minimal data for the abundance report."""
        np.random.seed(42)
        results = pd.DataFrame({
            'feature': ['Ala', 'Gly', 'Lys', 'Pro', 'Ser'],
            'log2FoldChange': [2.1, -1.5, 0.3, -2.8, 1.2],
            'padj': [0.001, 0.02, 0.8, 0.0001, 0.04],
            'baseMean': [500, 300, 150, 800, 250],
            'comparison': ['KO_vs_WT'] * 5,
        })

        count_matrix = pd.DataFrame(
            np.random.randint(50, 500, (4, 5)),
            index=['s1', 's2', 's3', 's4'],
            columns=['Ala', 'Gly', 'Lys', 'Pro', 'Ser'],
        )

        condition_map = {'s1': 'WT', 's2': 'WT', 's3': 'KO', 's4': 'KO'}

        return results, count_matrix, condition_map

    def test_report_generation(self, report_data, tmp_path):
        from trnaseq.qc.abundance_report import AbundanceReportGenerator
        results, count_matrix, condition_map = report_data

        gen = AbundanceReportGenerator(
            results_df=results,
            count_matrix=count_matrix,
            control_group='WT',
            level='aa',
            condition_map=condition_map,
        )
        out = tmp_path / 'abundance_report.html'
        result = gen.generate_html_report(out)
        assert Path(result).exists()
        content = Path(result).read_text()
        assert 'tRNA Abundance Dashboard' in content

    def test_all_panels_present(self, report_data, tmp_path):
        from trnaseq.qc.abundance_report import AbundanceReportGenerator
        results, count_matrix, condition_map = report_data

        gen = AbundanceReportGenerator(
            results_df=results,
            count_matrix=count_matrix,
            control_group='WT',
            level='aa',
            condition_map=condition_map,
        )
        out = tmp_path / 'abundance_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        assert 'Volcano Plot' in content
        assert 'MA Plot' in content
        assert 'Top Differentially Expressed' in content
        assert 'PCA' in content

    def test_empty_results(self, tmp_path):
        from trnaseq.qc.abundance_report import AbundanceReportGenerator
        gen = AbundanceReportGenerator(
            results_df=pd.DataFrame(),
            count_matrix=pd.DataFrame(),
            control_group='WT',
            level='aa',
        )
        out = tmp_path / 'abundance_report.html'
        result = gen.generate_html_report(out)
        assert Path(result).exists()

    def test_pca_with_few_samples(self, tmp_path):
        """PCA should be skipped if fewer than 3 samples."""
        from trnaseq.qc.abundance_report import AbundanceReportGenerator
        count_matrix = pd.DataFrame(
            {'Ala': [100, 200], 'Gly': [150, 300]},
            index=['s1', 's2'],
        )
        gen = AbundanceReportGenerator(
            results_df=pd.DataFrame(),
            count_matrix=count_matrix,
            control_group='WT',
            level='aa',
        )
        out = tmp_path / 'abundance_report.html'
        gen.generate_html_report(out)
        content = out.read_text()
        # PCA panel should not be present (need >= 3 samples)
        assert 'PCA' not in content


class TestAbundanceExport:
    """Test DifferentialAbundance.export_results (without DESeq2)."""

    @pytest.fixture
    def stats_csv(self, tmp_path):
        data = {
            'sample_name_unique': ['s1'] * 2 + ['s2'] * 2,
            'amino_acid': ['Ala', 'Gly'] * 2,
            'tRNA_annotation': ['tRNA-Ala', 'tRNA-Gly'] * 2,
            'codon': ['GCT', 'GGC'] * 2,
            'count': [100, 200, 120, 180],
        }
        csv_path = tmp_path / 'stats.csv'
        pd.DataFrame(data).to_csv(csv_path, index=False)
        return str(csv_path)

    @pytest.fixture
    def sample_df(self):
        return pd.DataFrame({
            'sample_name_unique': ['s1', 's2'],
            'sample_name': ['WT', 'KO'],
        })

    def test_count_matrix_export(self, stats_csv, sample_df, tmp_path):
        from trnaseq.abundance import DifferentialAbundance
        da = DifferentialAbundance(stats_csv, sample_df, level='aa')
        # Can't run DESeq2 without pydeseq2, but count matrix should work
        cm_path = tmp_path / 'count_matrix.csv'
        da.count_matrix.to_csv(cm_path)
        loaded = pd.read_csv(cm_path, index_col=0)
        assert loaded.shape == da.count_matrix.shape
