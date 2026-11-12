"""
Tests for per-sample stage 0c/1 stats files (trnaseq/io/sample_stats.py) and
their use by the pipeline runner.

SLURM array tasks run stages 0c/1 for one sample each and skip save_outputs,
so before these files existed the aggregation run never saw N_after_trim,
N_UMI_*, Mapping_percent etc. These tests cover: writing and merging the
files, the aggregation run's completeness check, array-mode row updates, the
stage 2 Mapping_percent denominator, the downsample_percentile guard, and the
shared sample_stats.xlsx switch on UMI_trim / SWIPE_align.

Run tests with:
    python -m pytest tests/test_sample_stats.py -v
"""

import json

import numpy as np
import pandas as pd
import pytest
import yaml

from src.alignment import SWIPE_align
from src.read_processing import UMI_trim
from trnaseq.cli.run_pipeline import CheckStatus, PreflightReport, PreprocessingPipeline
from trnaseq.io.sample_stats import (
    SAMPLE_STATS_SUBDIR,
    assign_sample_values,
    format_array_indices,
    merge_sample_stats,
    read_sample_stats,
    recompute_mapping_percent,
    sample_stats_path,
    write_sample_stats,
)


SAMPLES = ['S0', 'S1', 'S2']


def _stage0b_df():
    """sample_df as stage 0b leaves it: no 0c/1 columns."""
    return pd.DataFrame({'sample_name_unique': SAMPLES,
                         'species': ['ecoli'] * 3,
                         'N_total': [1000, 2000, 3000]})


def _write_0c(project, sample, n_trim, n_ds=None):
    n_ds = n_trim if n_ds is None else n_ds
    df = pd.DataFrame([{
        'sample_name_unique': sample, 'N_total': 9999,  # not a 0c column: must not be written
        'N_after_trim': n_trim, 'N_UMI_observed': 50, 'N_UMI_expected': 60.5,
        'N_after_downsample': n_ds, 'percent_seqs_after_UMI_trim': 80.0,
        'percent_UMI_obs-vs-exp': 82.6,
    }])
    return write_sample_stats(project, df, '0c')


def _write_1(project, sample, n_mapped, mapping_pct):
    df = pd.DataFrame([{
        'sample_name_unique': sample, 'N_mapped': n_mapped,
        'percent_single_annotation': 90.0, 'percent_multiple_annotation': 10.0,
        'percent_multiple_codons': 1.0, 'Mapping_percent': mapping_pct,
    }])
    return write_sample_stats(project, df, '1')


# ---------------------------------------------------------------------------
# Writing / reading files
# ---------------------------------------------------------------------------

def test_write_only_stage_columns_and_no_temp_files(tmp_path):
    (path,) = _write_0c(tmp_path, 'S0', 700)

    assert path == sample_stats_path(tmp_path, 'S0', '0c')
    assert path.parent == tmp_path / SAMPLE_STATS_SUBDIR
    payload = read_sample_stats(path, 'S0', '0c')
    assert payload['values']['N_after_trim'] == 700
    assert 'N_total' not in payload['values']
    assert set(payload['provenance']) >= {'written_at', 'host', 'slurm_array_task_id'}
    assert [p.name for p in path.parent.iterdir()] == ['S0.stage0c.json']


def test_numpy_and_nan_values_round_trip(tmp_path):
    df = pd.DataFrame({'sample_name_unique': ['S0'],
                       'N_after_trim': np.array([7], dtype=np.int64),
                       'N_after_downsample': [np.nan]})
    (path,) = write_sample_stats(tmp_path, df, '0c')

    raw = json.loads(path.read_text())
    assert raw['values'] == {'N_after_trim': 7, 'N_after_downsample': None}

    merged, _ = merge_sample_stats(_stage0b_df(), tmp_path, stages=['0c'])
    assert merged.loc[0, 'N_after_trim'] == 7
    assert np.isnan(merged.loc[0, 'N_after_downsample'])


def test_file_for_wrong_sample_is_rejected(tmp_path):
    (path,) = _write_0c(tmp_path, 'S0', 700)
    path.rename(sample_stats_path(tmp_path, 'S1', '0c'))

    with pytest.raises(ValueError, match="holds sample 'S0'"):
        merge_sample_stats(_stage0b_df(), tmp_path, stages=['0c'])


# ---------------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------------

def test_merge_two_array_tasks_keeps_row_order_and_reports_gap(tmp_path):
    _write_0c(tmp_path, 'S2', 2500)
    _write_0c(tmp_path, 'S0', 700)

    merged, status = merge_sample_stats(_stage0b_df(), tmp_path)

    assert merged['sample_name_unique'].tolist() == SAMPLES
    assert merged['N_total'].tolist() == [1000, 2000, 3000]
    assert merged.loc[0, 'N_after_trim'] == 700
    assert np.isnan(merged.loc[1, 'N_after_trim'])
    assert merged.loc[2, 'N_after_trim'] == 2500
    assert status['0c'].n_files == 2
    assert status['0c'].missing == [(1, 'S1')]
    assert len(status['1'].missing) == 3


def test_file_values_override_stale_sample_df(tmp_path):
    df = _stage0b_df()
    df['N_after_trim'] = [111, 222, 333]
    df['N_mapped'] = [1, 2, 3]
    _write_0c(tmp_path, 'S1', 1800)

    merged, status = merge_sample_stats(df, tmp_path)

    assert merged['N_after_trim'].tolist() == [111, 1800, 333]
    # Values already in sample_df count as coverage, so nothing is missing:
    assert status['0c'].complete and status['1'].complete


def test_merge_restricted_to_samples(tmp_path):
    for s in SAMPLES:
        _write_0c(tmp_path, s, 500)

    merged, status = merge_sample_stats(_stage0b_df(), tmp_path, samples=['S1'])

    assert merged['N_after_trim'].isna().tolist() == [True, False, True]
    assert status['0c'].n_samples == 1 and status['0c'].complete


def test_assign_sample_values_adds_and_overwrites():
    df = pd.DataFrame({'sample_name_unique': SAMPLES, 'a': [1.0, 2.0, 3.0]})
    values = pd.DataFrame({'a': [20.0], 'b': [5.0]}, index=['S1'])

    out = assign_sample_values(df, values)

    assert out['a'].tolist() == [1.0, 20.0, 3.0]
    assert out['b'].isna().tolist() == [True, False, True]
    assert df['a'].tolist() == [1.0, 2.0, 3.0]  # input untouched


@pytest.mark.parametrize('idx,expected', [
    ([], ''), ([3], '3'), ([0, 1, 2, 5, 7, 8], '0-2,5,7-8'), ([4, 2, 3, 3], '2-4'),
])
def test_format_array_indices(idx, expected):
    assert format_array_indices(idx) == expected


# ---------------------------------------------------------------------------
# Stage 2 Mapping_percent denominator
# ---------------------------------------------------------------------------

def test_recompute_mapping_percent_uses_downsampled_count():
    df = pd.DataFrame({
        'sample_name_unique': ['ds', 'no_ds_col', 'zero', 'legacy'],
        'N_total_aligned': [400, 300, 0, 50],
        'N_after_trim': [1000, 600, 0, np.nan],
        'N_after_downsample': [500, np.nan, 0, np.nan],
        'Mapping_percent': [1.0, 1.0, 1.0, 77.0],
    })

    out = recompute_mapping_percent(df)

    assert out['Mapping_percent'].tolist() == pytest.approx([80.0, 50.0, 0.0, 77.0])


def test_recompute_mapping_percent_without_denominator_columns():
    df = pd.DataFrame({'sample_name_unique': ['S0'], 'N_total_aligned': [5],
                       'Mapping_percent': [42.0]})
    assert recompute_mapping_percent(df)['Mapping_percent'].tolist() == [42.0]


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

@pytest.fixture
def make_pipeline(tmp_path):
    def _make(sample_index=None, **config):
        cfg = {'sample_list': 'sample_list.xlsx', 'run_abundance_analysis': False}
        cfg.update(config)
        cfg_path = tmp_path / 'config.yaml'
        cfg_path.write_text(yaml.safe_dump(cfg))
        pipe = PreprocessingPipeline(cfg_path, tmp_path, n_jobs=1, sample_index=sample_index)
        pipe.sample_df = _stage0b_df()
        return pipe
    return _make


def test_gather_run_merges_files_into_sample_df_xlsx(tmp_path, make_pipeline):
    """End to end through run(): load sample_df.xlsx, merge, save."""
    _stage0b_df().to_excel(tmp_path / 'sample_df.xlsx', index=False)
    for i, s in enumerate(SAMPLES):
        _write_0c(tmp_path, s, 1000 * (i + 1))
        _write_1(tmp_path, s, 100 * (i + 1), 10.0 * (i + 1))

    make_pipeline().run(stages={'7'})

    saved = pd.read_excel(tmp_path / 'sample_df.xlsx')
    assert saved['sample_name_unique'].tolist() == SAMPLES
    assert saved['N_after_trim'].tolist() == [1000, 2000, 3000]
    assert saved['N_mapped'].tolist() == [100, 200, 300]
    assert saved['Mapping_percent'].tolist() == [10.0, 20.0, 30.0]
    assert not list(tmp_path.glob('.sample_df.*'))


def test_gather_run_fails_fast_on_missing_samples(tmp_path, make_pipeline):
    _stage0b_df().to_excel(tmp_path / 'sample_df.xlsx', index=False)
    for s in ('S0', 'S2'):
        _write_0c(tmp_path, s, 1000)
        _write_1(tmp_path, s, 100, 10.0)

    with pytest.raises(RuntimeError) as exc:
        make_pipeline().run(stages={'7'})

    msg = str(exc.value)
    assert 'stage 0c: 1/3 samples missing (array indices 1): S1' in msg
    assert 'stage 1: 1/3 samples missing' in msg
    assert 'sbatch --array=1 ' in msg


def test_allow_incomplete_sample_stats_warns_and_saves(make_pipeline, tmp_path):
    _write_0c(tmp_path, 'S0', 1000)
    pipe = make_pipeline(allow_incomplete_sample_stats=True)
    pipe.sample_df['N_mapped'] = [1, 2, 3]

    pipe._merge_sample_stats({'2', '3'})

    assert (tmp_path / 'sample_df.xlsx').exists()
    assert 'Per-sample stats are incomplete' in (tmp_path / 'logs' / 'pipeline.log').read_text()


def test_legacy_project_without_any_files_only_warns(make_pipeline, tmp_path):
    pipe = make_pipeline()

    pipe._merge_sample_stats({'2', '3', '5'})

    assert 'N_after_trim' not in pipe.sample_df.columns
    assert not (tmp_path / 'sample_df.xlsx').exists()
    log = (tmp_path / 'logs' / 'pipeline.log').read_text()
    assert 'No stage 0c stats for any sample' in log
    assert 'No stage 1 stats for any sample' in log


def test_no_completeness_check_for_stages_this_run_produces(make_pipeline, tmp_path):
    _write_0c(tmp_path, 'S0', 1000)
    pipe = make_pipeline()

    pipe._merge_sample_stats({'0c', '1', '2'})   # would raise if 0c/1 were checked
    pipe._merge_sample_stats({'0c'})             # no downstream stage: no check


def test_stage_0b_runs_skip_merge(make_pipeline, tmp_path):
    for s in SAMPLES:
        _write_0c(tmp_path, s, 1000)
    pipe = make_pipeline()

    pipe._merge_sample_stats({'0a', '0b', '0c', '1', '2'})

    assert 'N_after_trim' not in pipe.sample_df.columns


def test_array_task_merges_only_its_own_sample(make_pipeline, tmp_path):
    for s in SAMPLES:
        _write_0c(tmp_path, s, 1000)
    pipe = make_pipeline(sample_index=1)

    pipe._merge_sample_stats({'1'})   # S0/S2 lack stage 1 too: must not raise

    assert pipe.sample_df['N_after_trim'].isna().tolist() == [True, False, True]
    assert not (tmp_path / 'sample_df.xlsx').exists()


def test_record_stage_stats_array_mode_overwrites_stale_row(make_pipeline, tmp_path):
    pipe = make_pipeline(sample_index=1)
    pipe.sample_df['N_after_trim'] = [111, 222, 333]
    result_df = pd.DataFrame([{'sample_name_unique': 'S1', 'species': 'ecoli',
                               'N_total': 2000, 'N_after_trim': 1800,
                               'N_after_downsample': 1800}])

    pipe._record_stage_stats('0c', result_df)

    assert pipe.sample_df['N_after_trim'].tolist() == [111, 1800, 333]
    assert pipe.sample_df.loc[1, 'N_after_downsample'] == 1800
    assert [p.name for p in (tmp_path / SAMPLE_STATS_SUBDIR).iterdir()] == ['S1.stage0c.json']


# ---------------------------------------------------------------------------
# downsample_percentile under --sample-index
# ---------------------------------------------------------------------------

def test_stage_0c_rejects_downsample_percentile_in_array_mode(make_pipeline):
    pipe = make_pipeline(sample_index=0, downsample_percentile=50)
    with pytest.raises(ValueError, match='downsample_percentile is not supported'):
        pipe.stage_0c_umi_trim()


@pytest.mark.parametrize('sample_index,expected', [(0, CheckStatus.FAIL), (None, CheckStatus.WARN)])
def test_preflight_flags_downsample_percentile(make_pipeline, sample_index, expected):
    pipe = make_pipeline(sample_index=sample_index, downsample_percentile=50)
    report = PreflightReport()

    pipe._preflight_config(report, {'0c', '1'})

    (res,) = [r for r in report.results if r.name == 'downsample_percentile']
    assert res.status == expected


def test_preflight_silent_without_downsample_percentile(make_pipeline):
    report = PreflightReport()
    make_pipeline(sample_index=0)._preflight_config(report, {'0c', '1'})
    assert not [r for r in report.results if r.name == 'downsample_percentile']


# ---------------------------------------------------------------------------
# Shared sample_stats.xlsx switch
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('flag', [True, False])
def test_umi_trim_write_stats_xlsx(tmp_path, flag):
    trim = UMI_trim.__new__(UMI_trim)
    trim.UMI_dir_abs = str(tmp_path)
    trim.write_stats_xlsx = flag
    trim.sample_df = pd.DataFrame({'sample_name_unique': ['S0'], 'N_total': [10]})

    trim._collect_stats([['S0', 8, 5, 6.0, 8]])

    assert (tmp_path / 'sample_stats.xlsx').exists() is flag
    assert trim.sample_df.loc[0, 'N_after_downsample'] == 8


@pytest.mark.parametrize('flag', [True, False])
def test_swipe_align_write_stats_xlsx(tmp_path, monkeypatch, flag):
    monkeypatch.chdir(tmp_path)
    aligner = SWIPE_align.__new__(SWIPE_align)
    aligner.write_stats_xlsx = flag
    aligner.sample_df = pd.DataFrame({'sample_name_unique': ['S0']})

    aligner._write_stats([['S0', 5, 100.0, 0.0, 0.0, 62.5]])

    assert (tmp_path / 'sample_stats.xlsx').exists() is flag
    assert aligner.sample_df.loc[0, 'Mapping_percent'] == 62.5
