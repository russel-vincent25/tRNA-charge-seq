"""
Regression tests for :meth:`SWIPE_align._collect_stats` mapping percentage.

``N_after_downsample`` is created in memory by stage 0c (UMI trimming) and
is absent from a sample_df loaded from disk when a run starts at stage 1
(e.g. SLURM array tasks, which skip ``save_outputs``). ``_collect_stats``
must then fall back to ``Ninput``, the read count of the UMI-trimmed file
that SWIPE aligned, instead of raising ``KeyError``.

Run tests with:
    python -m pytest tests/test_alignment_stats.py -v
"""

import bz2
import json

import numpy as np
import pandas as pd
import pytest

from src.alignment import SWIPE_align


SAMPLE = 'S1'
N_READS = 8
N_MAPPED = 5


def _make_aligner(umi_dir):
    """Bare SWIPE_align with only the attributes _collect_stats reads."""
    aligner = SWIPE_align.__new__(SWIPE_align)
    aligner.verbose = False
    aligner.stream = False
    aligner.common_seqs_fnam = None
    aligner.from_UMIdir = True
    aligner.UMI_dir_abs = str(umi_dir)
    return aligner


@pytest.fixture
def align_dir(tmp_path, monkeypatch):
    """Alignment dir (cwd) with a SWalign JSON plus a UMI-trimmed FASTQ."""
    umi_dir = tmp_path / 'UMI_trimmed'
    umi_dir.mkdir()
    with bz2.open(umi_dir / f'{SAMPLE}_UMI-trimmed.fastq.bz2', 'wt') as fh:
        for i in range(N_READS):
            fh.write(f'@read{i}\nACGTACGT\n+\nIIIIIIII\n')

    align_dir = tmp_path / 'SWIPE_alignment'
    align_dir.mkdir()
    swres = {}
    for i in range(N_READS):
        if i < N_MAPPED:
            swres[f'read{i}'] = {'aligned': True,
                                 'name': 'Ecoli_tRNA-Ala-AGC-1-1',
                                 'one_codon': True}
        else:
            swres[f'read{i}'] = {'aligned': False}
    with bz2.open(align_dir / f'{SAMPLE}_SWalign.json.bz2', 'wt') as fh:
        json.dump(swres, fh)

    monkeypatch.chdir(align_dir)
    return umi_dir


def _row(sample_df):
    return next(sample_df.iterrows())


def test_mapping_percent_without_n_after_downsample(align_dir):
    sample_df = pd.DataFrame({'sample_name_unique': [SAMPLE]})
    assert 'N_after_downsample' not in sample_df.columns

    res = _make_aligner(align_dir)._collect_stats(*_row(sample_df))

    assert res[0] == SAMPLE
    assert res[1] == N_MAPPED
    assert res[-1] == pytest.approx(N_MAPPED / N_READS * 100)


def test_mapping_percent_with_nan_n_after_downsample(align_dir):
    sample_df = pd.DataFrame({'sample_name_unique': [SAMPLE],
                              'N_after_downsample': [np.nan]})

    res = _make_aligner(align_dir)._collect_stats(*_row(sample_df))

    assert res[-1] == pytest.approx(N_MAPPED / N_READS * 100)


def test_mapping_percent_uses_n_after_downsample_when_present(align_dir):
    sample_df = pd.DataFrame({'sample_name_unique': [SAMPLE],
                              'N_after_downsample': [10]})

    res = _make_aligner(align_dir)._collect_stats(*_row(sample_df))

    assert res[-1] == pytest.approx(N_MAPPED / 10 * 100)


def test_mapping_percent_zero_guard(align_dir):
    sample_df = pd.DataFrame({'sample_name_unique': [SAMPLE],
                              'N_after_downsample': [0]})

    res = _make_aligner(align_dir)._collect_stats(*_row(sample_df))

    assert res[-1] == 0
