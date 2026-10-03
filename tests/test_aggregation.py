"""Tests for site-level replicate aggregation (CHANNEL_FIX_SPEC change 5)."""

import numpy as np
import pytest
import pandas as pd

from trnaseq.modifications.modification_caller import ReplicateAggregator


def _call(mod='m1A', source='known_modomics', mm=0.3, gap=0.05, stop=20.0,
          pattern='A->T', channel='mismatch', fired='mismatch', pval=1e-8, **extra):
    row = {'trna_name': 'tRNA-Ala-AGC-1-1', 'position': 58, 'modified_position': 58,
           'modification': mod, 'source': source, 'identity_support': 'modomics',
           'mismatch_rate': mm, 'gap_rate': gap, 'rt_stop_pct': stop,
           'dominant_pattern': pattern, 'dominant_channel': channel,
           'channels_fired': fired, 'fold_change': 300.0, 'coverage': 1000,
           'confidence': 0.8, 'pvalue': pval, 'rt_enzyme': 'Maxima'}
    row.update(extra)
    return pd.DataFrame([row])


GROUPS = {'WT': ['S1', 'S2', 'S3']}


def test_differing_labels_stay_one_site():
    per_sample = {
        'S1': _call(),
        'S2': _call(),
        'S3': _call(mod='novel_candidate', source='novel_candidate',
                    channel='deletion', fired='deletion', pattern=''),
    }
    out = ReplicateAggregator(min_replicates=3).aggregate(per_sample, GROUPS)
    assert len(out) == 1
    row = out.iloc[0]
    assert row['n_replicates_detected'] == 3 and row['consensus_call']
    assert row['modification'] == 'm1A'
    assert row['labels'] == 'm1A:2; novel_candidate:1'
    assert row['channels_fired'] == 'mismatch:2; deletion:1'


def test_per_channel_rates_and_pattern_carried_through():
    per_sample = {'S1': _call(mm=0.2, gap=0.1, stop=10.0),
                  'S2': _call(mm=0.4, gap=0.3, stop=30.0, pattern='A->G'),
                  'S3': _call(mm=0.3, gap=0.2, stop=20.0)}
    row = ReplicateAggregator().aggregate(per_sample, GROUPS).iloc[0]
    assert row['mean_mismatch_rate'] == pytest.approx(0.3)
    assert row['mean_gap_rate'] == pytest.approx(0.2)
    assert row['mean_rt_stop_pct'] == pytest.approx(20.0)
    assert row['dominant_pattern'] == 'A->T'
    assert row['dominant_channel'] == 'mismatch'
    assert row['modified_position'] == 58 and row['rt_enzyme'] == 'Maxima'


def test_novel_only_when_no_replicate_resolves_identity():
    novel = dict(mod='novel_candidate', source='novel_candidate')
    out = ReplicateAggregator().aggregate(
        {'S1': _call(**novel), 'S2': _call(**novel)}, GROUPS)
    assert out.iloc[0]['modification'] == 'novel_candidate'
    assert out.iloc[0]['source'] == 'novel_candidate'


def test_one_pvalue_per_replicate():
    # Legacy tables can hold several rows per site and sample
    per_sample = {'S1': pd.concat([_call(pval=1e-3), _call(mod='m6A', pval=1e-9)])}
    row = ReplicateAggregator(min_replicates=1).aggregate(per_sample, GROUPS).iloc[0]
    assert row['n_replicates_detected'] == 1
    assert row['fisher_combined_pvalue'] == 1e-9


def test_missing_channel_columns_are_nan():
    legacy = pd.DataFrame([{'trna_name': 't', 'position': 1, 'modification': 'm1A',
                            'mismatch_rate': 0.2, 'coverage': 100,
                            'confidence': 0.5, 'pvalue': 1e-5}])
    row = ReplicateAggregator(min_replicates=1).aggregate({'S1': legacy}, GROUPS).iloc[0]
    assert np.isnan(row['mean_gap_rate']) and np.isnan(row['mean_rt_stop_pct'])
    assert row['dominant_pattern'] == ''
