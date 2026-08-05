"""
Shared pytest fixtures for the trnaseq test suite.

The :func:`synth_run` fixture is the workhorse: a single namespace of
small (≤50-sample) synthetic DataFrames + a :class:`ReportContext`,
sufficient to exercise every panel across every QC report without
touching real pipeline outputs. Per PRP §4.6 individual tests should
finish in <30 s of wall time, so the data is intentionally tiny.
"""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_AAS = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His',
        'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp',
        'Tyr', 'Val']

_ANTICODONS = ['AGC', 'CGT', 'ATT', 'GTC', 'GCA', 'CTG', 'TTC', 'GCC',
               'GTG', 'GAT', 'CAG', 'CTT', 'CAT', 'GAA', 'CGG', 'TGA',
               'AGT', 'CCA', 'GTA', 'AAC']


def _make_trna_names(n: int) -> list[str]:
    """Return E. coli-style tRNA annotation names."""
    return [f'Ecoli-{_AAS[i % len(_AAS)]}-{_ANTICODONS[i % len(_ANTICODONS)]}-1-{(i // 10) + 1}'
            for i in range(n)]


# ---------------------------------------------------------------------------
# Bundle
# ---------------------------------------------------------------------------

@dataclass
class SynthRun:
    """One bag holding every DataFrame the generator tests need.

    Field names mirror the kwargs the generator constructors expect, so a
    test can splat where convenient::

        report = QCReportGenerator(
            project_dir=tmp_path,
            sample_df=synth_run.sample_df,
            inp_file_df=synth_run.inp_file_df,
            ...
        )
    """
    project_dir: Path
    sample_df: pd.DataFrame
    inp_file_df: pd.DataFrame
    charge_summary_df: pd.DataFrame
    charge_df_transcript: pd.DataFrame
    charge_df_aa: pd.DataFrame
    fragment_counts_df: pd.DataFrame
    rt_dropoff_df: pd.DataFrame
    fragment_lengths_df: pd.DataFrame
    fragment_summary_df: pd.DataFrame
    stats_df: pd.DataFrame
    abundance_results_df: pd.DataFrame
    abundance_count_matrix: pd.DataFrame
    abundance_condition_map: dict
    context: object  # ReportContext, kept untyped to avoid import-time cost


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope='session')
def n_samples() -> int:
    """50 samples — comfortably above the legacy 30-pair cutoff so panels
    exercise their large-N branches."""
    return 50


@pytest.fixture(scope='session')
def synth_run(tmp_path_factory, n_samples) -> SynthRun:  # noqa: D401
    """Construct a synthetic pipeline output suitable for all 5 generators."""
    from trnaseq.qc._common import ReportContext

    rng = np.random.default_rng(42)
    project_dir = tmp_path_factory.mktemp('synth_project')

    samples = [f'S{i:03d}' for i in range(n_samples)]
    # 10 replicate groups → 5 reps each
    conditions = [f'cond_{i // 5}' for i in range(n_samples)]

    sample_df = pd.DataFrame({
        'sample_name_unique': samples,
        'sample_name': conditions,
        'species': ['E. coli'] * n_samples,
        'fastq_mate1_filename': [f'lane01_{s}.fastq.gz' for s in samples],
        'N_total': rng.integers(50_000, 200_000, n_samples),
        'N_after_trim': rng.integers(40_000, 180_000, n_samples),
        'N_UMI_observed': rng.integers(8_000, 30_000, n_samples),
        'N_UMI_expected': rng.integers(10_000, 35_000, n_samples),
        'Mapping_percent': rng.uniform(35, 95, n_samples),
        'percent_seqs_after_UMI_trim': rng.uniform(60, 99, n_samples),
        'percent_single_annotation': rng.uniform(50, 95, n_samples),
        'N_full_length': rng.integers(20_000, 80_000, n_samples),
        'N_rt_dropoff': rng.integers(2_000, 20_000, n_samples),
        'N_5p_fragment': rng.integers(500, 5_000, n_samples),
        'N_degraded': rng.integers(100, 2_000, n_samples),
        'N_total_aligned': rng.integers(40_000, 150_000, n_samples),
    })

    inp_file_df = pd.DataFrame({
        'fastq_mate1_filename': [f'lane01_{s}.fastq.gz' for s in samples],
        'N_pairs': rng.integers(100_000, 400_000, n_samples),
        'N_merged': rng.integers(90_000, 380_000, n_samples),
        'percent_successfully_merged': rng.uniform(70, 99, n_samples),
        'percent_BC-mapped': rng.uniform(50, 95, n_samples),
    })

    # --- charge dataframes -------------------------------------------------
    n_trnas = 40
    trna_names = _make_trna_names(n_trnas)
    charge_summary = pd.DataFrame({
        'sample_name_unique': samples,
        'charge_canonical_mean': rng.uniform(25, 90, n_samples),
        'tRNA_source': ['host'] * n_samples,
    })

    long_rows = []
    for s in samples:
        for t in trna_names:
            aa = t.split('-')[1]
            long_rows.append({
                'sample_name_unique': s,
                'tRNA_annotation': t,
                'amino_acid': aa,
                'charge_canonical': float(rng.uniform(20, 95)),
                'count': int(rng.integers(50, 5_000)),
                'RPM': float(rng.uniform(10, 5_000)),
                'tRNA_source': 'host',
            })
    charge_df_transcript = pd.DataFrame(long_rows)

    aa_rows = []
    for s in samples:
        for aa in sorted(set(t.split('-')[1] for t in trna_names)):
            aa_rows.append({
                'sample_name_unique': s,
                'amino_acid': aa,
                'charge_canonical': float(rng.uniform(20, 95)),
                'tRNA_source': 'host',
            })
    charge_df_aa = pd.DataFrame(aa_rows)

    # --- fragment dataframes ----------------------------------------------
    frag_rows = []
    for s in samples:
        for t in trna_names:
            frag_rows.append({
                'sample_name_unique': s,
                'tRNA_annotation': t,
                'integrity_score': float(rng.uniform(0.5, 1.0)),
                'N_full_length': int(rng.integers(20, 500)),
                'N_rt_dropoff': int(rng.integers(5, 100)),
                'N_5p_tRF': int(rng.integers(1, 30)),
                'N_degraded': int(rng.integers(0, 10)),
            })
    fragment_counts_df = pd.DataFrame(frag_rows)

    rt_drop_rows = []
    for s in samples:
        for pos in range(0, 80, 5):
            rt_drop_rows.append({
                'sample_name_unique': s,
                'position': pos,
                'rt_stop_fraction': float(rng.uniform(0, 0.1)),
            })
    rt_dropoff_df = pd.DataFrame(rt_drop_rows)

    frag_len_rows = []
    for s in samples:
        for length in range(20, 90, 10):
            frag_len_rows.append({
                'sample_name_unique': s,
                'fragment_length': length,
                'count': int(rng.integers(100, 5_000)),
            })
    fragment_lengths_df = pd.DataFrame(frag_len_rows)

    fragment_summary_df = sample_df[[
        'sample_name_unique', 'N_full_length',
        'N_rt_dropoff', 'N_5p_fragment', 'N_degraded',
    ]].copy().rename(columns={'N_5p_fragment': 'N_5p_tRF'})

    # --- stats_df (long-form, used by PCA/replicate panels) ---------------
    stats_rows = []
    for s in samples:
        for t in trna_names:
            stats_rows.append({
                'sample_name_unique': s,
                'tRNA_annotation': t,
                'count': int(rng.integers(50, 5_000)),
            })
    stats_df = pd.DataFrame(stats_rows)

    # --- abundance dataframes ---------------------------------------------
    cond_pairs = [('cond_0', f'cond_{i}') for i in range(1, 5)]
    abundance_rows = []
    for ctrl, treat in cond_pairs:
        for t in trna_names:
            abundance_rows.append({
                'feature': t,
                'baseMean': float(rng.uniform(50, 5000)),
                'log2FoldChange': float(rng.normal(0, 1.5)),
                'padj': float(rng.uniform(1e-6, 1)),
                'comparison': f'{treat}_vs_{ctrl}',
            })
    abundance_results_df = pd.DataFrame(abundance_rows)

    # Counts matrix: samples x features
    counts_mat = rng.integers(50, 5000, size=(n_samples, n_trnas))
    abundance_count_matrix = pd.DataFrame(counts_mat, index=samples, columns=trna_names)

    abundance_condition_map = dict(zip(samples, conditions))

    ctx = ReportContext(
        project_dir=project_dir,
        pipeline_version='0.0.0-test',
        config_path=project_dir / 'config.yaml',
        config_sha256='a' * 64,
        sample_sheet_path=project_dir / 'sample_sheet.xlsx',
        sample_sheet_sha256='b' * 64,
        reference_db='ecoli_test (20 sequences)',
        generated_at=_dt.datetime(2026, 5, 27, 14, 32, 11, tzinfo=_dt.timezone.utc),
        runtime_seconds=42.0,
        host='test-host',
        command='pytest tests/qc/',
    )

    return SynthRun(
        project_dir=project_dir,
        sample_df=sample_df,
        inp_file_df=inp_file_df,
        charge_summary_df=charge_summary,
        charge_df_transcript=charge_df_transcript,
        charge_df_aa=charge_df_aa,
        fragment_counts_df=fragment_counts_df,
        rt_dropoff_df=rt_dropoff_df,
        fragment_lengths_df=fragment_lengths_df,
        fragment_summary_df=fragment_summary_df,
        stats_df=stats_df,
        abundance_results_df=abundance_results_df,
        abundance_count_matrix=abundance_count_matrix,
        abundance_condition_map=abundance_condition_map,
        context=ctx,
    )
