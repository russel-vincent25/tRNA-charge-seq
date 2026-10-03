"""Tests for the novel-site summary and its report panel."""

import numpy as np
import pandas as pd

from trnaseq.modifications.novel_sites import summarize_novel_sites
from trnaseq.qc.modification_report import ModificationReportGenerator


def _call(position, source='novel_candidate', modification='novel_candidate',
          channels='deletion', dominant='deletion', pattern='', modified=None,
          enzyme='Maxima', **rates):
    row = {'trna_name': 'tRNA-Pro-TGG-1-1', 'position': position,
           'modified_position': position if modified is None else modified,
           'source': source, 'modification': modification,
           'channels_fired': channels, 'dominant_channel': dominant,
           'dominant_pattern': pattern, 'candidates': 'm7G;m1G',
           'mismatch_rate': 0.0, 'gap_rate': 0.0, 'rt_stop_pct': 0.0,
           'coverage': 1000, 'rt_enzyme': enzyme}
    row.update(rates)
    return row


def _per_sample():
    return {
        'Maxima-55-R1': pd.DataFrame([_call(36, gap_rate=0.12),
                                      _call(38, source='known_modomics', modification='m1G')]),
        'Maxima-55-R2': pd.DataFrame([_call(36, gap_rate=0.10)]),
        'SSIV-55-R1': pd.DataFrame([_call(36, channels='mismatch+deletion', pattern='G->C',
                                          gap_rate=0.2, mismatch_rate=0.15, enzyme='SSIV'),
                                    _call(60, channels='rt_stop', dominant='rt_stop',
                                          modified=59, rt_stop_pct=40.0, enzyme='SSIV')]),
    }


GROUPS = {'Maxima-55': ['Maxima-55-R1', 'Maxima-55-R2'], 'SSIV-55': ['SSIV-55-R1']}
KNOWN = {'tRNA-Pro-TGG-1-1': pd.DataFrame({'linear_position': [38],
                                           'modification_short_name': ['m1G']})}
REF = {'tRNA-Pro-TGG-1-1': {'seq': 'G' * 76}}


def test_one_row_per_novel_site():
    ns = summarize_novel_sites(_per_sample(), GROUPS, KNOWN, REF)
    assert list(ns['position']) == [36, 60]          # most observed first
    row = ns.iloc[0]
    assert row['n_observations'] == 3 and row['n_samples'] == 3
    assert row['n_conditions'] == 2
    assert row['rt_enzymes'] == 'Maxima;SSIV'
    assert row['channels'] == 'deletion:2; mismatch+deletion:1'
    assert row['substitutions'] == 'G->C:1'
    assert row['median_gap_rate'] == 0.12
    assert row['ref_nt'] == 'G'


def test_nearest_known_offset_from_implicated_base():
    ns = summarize_novel_sites(_per_sample(), GROUPS, KNOWN, REF).set_index('position')
    assert ns.loc[36, 'nearest_known'] == 'm1G'
    assert ns.loc[36, 'nearest_known_offset'] == 2
    assert ns.loc[60, 'modified_position'] == 59
    assert ns.loc[60, 'nearest_known'] == '' and np.isnan(ns.loc[60, 'nearest_known_offset'])


def test_site_labelled_in_other_samples_is_reported():
    per = _per_sample()
    per['SSIV-55-R1'] = pd.concat([per['SSIV-55-R1'], pd.DataFrame(
        [_call(36, source='known_modomics', modification='m7G')])])
    ns = summarize_novel_sites(per, GROUPS, KNOWN, REF).set_index('position')
    assert ns.loc[36, 'labelled_elsewhere'] == 'm7G:1'


def test_no_novel_calls_gives_empty_table():
    per = {'s': pd.DataFrame([_call(38, source='known_modomics', modification='m1G')])}
    assert summarize_novel_sites(per).empty
    assert summarize_novel_sites({}).empty


def test_report_panel(tmp_path):
    per = _per_sample()
    ns = summarize_novel_sites(per, GROUPS, KNOWN, REF)
    out = ModificationReportGenerator(per, replicate_groups=GROUPS,
                                      novel_sites=ns).generate_html_report(
        tmp_path / 'modification_report.html')
    html = open(out).read()
    assert 'Novel Sites' in html
    assert '1 lie within 3 nt of a known modification' in html
    assert (tmp_path / 'data' / 'novel_sites.csv').exists()
