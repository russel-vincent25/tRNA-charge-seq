"""
Tests for channel-aware site calling (CHANNEL_FIX_SPEC change 4):
- Detection is OR across channels and invariant to which channel carries signal
- Per-channel beta-binomial backgrounds are calibrated on overdispersed nulls
- Channel priors change confidence and flags, never which sites are called
- One row per site with a ranked candidate list
- The analyzer's mismatch rate excludes deletions
"""

import numpy as np
import pandas as pd
import pytest

from trnaseq.modifications.modification_caller import (
    CHANNELS,
    ChannelBackground,
    ModificationCaller,
    _bh_qvalues,
    estimate_channel_backgrounds,
)
from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer

BACKGROUNDS = {
    'mismatch': ChannelBackground(mean=8e-4, rho=1e-3),
    'deletion': ChannelBackground(mean=1e-4, rho=5e-3),
    'rt_stop': ChannelBackground(mean=1e-3, rho=1e-2),
}


def _signatures(sites, n_null=40, coverage=2000):
    """Signature table: *n_null* clean positions plus {position: overrides}."""
    rows = []
    for pos in range(2, 2 + n_null):
        rows.append({'position': pos, 'coverage': coverage, 'correct_nt': 'G',
                     'mismatch_rate': 0.0008, 'gap_rate': 0.0001,
                     'rt_stop_pct': 0.1})
    for pos, over in sites.items():
        row = {'position': pos, 'coverage': coverage, 'correct_nt': 'G',
               'mismatch_rate': 0.0008, 'gap_rate': 0.0001, 'rt_stop_pct': 0.1}
        row.update(over)
        rows.append(row)
    return pd.DataFrame(rows)


def _caller(**kwargs):
    kwargs.setdefault('channel_backgrounds', BACKGROUNDS)
    return ModificationCaller(organism='Escherichia coli', **kwargs)


class TestChannelInvariance:

    @pytest.mark.parametrize('channel,override', [
        ('mismatch', {'mismatch_rate': 0.30}),
        ('deletion', {'gap_rate': 0.30}),
        ('rt_stop', {'rt_stop_pct': 40.0}),
    ])
    def test_signal_in_any_single_channel_is_called(self, channel, override):
        calls = _caller().call_all('tRNA-x', _signatures({100: override}))
        assert list(calls['position']) == [100]
        assert calls.iloc[0]['channels_fired'] == channel
        assert calls.iloc[0]['dominant_channel'] == channel

    def test_multiple_channels_recorded(self):
        calls = _caller().call_all(
            'tRNA-x', _signatures({100: {'gap_rate': 0.3, 'rt_stop_pct': 40.0}}))
        assert calls.iloc[0]['channels_fired'] == 'deletion+rt_stop'
        assert calls.iloc[0]['n_channels_fired'] == 2

    def test_below_threshold_not_called_even_if_significant(self):
        # 8% deletion is overwhelmingly significant but under the 10% threshold
        calls = _caller().call_all('tRNA-x', _signatures({100: {'gap_rate': 0.08}}))
        assert calls.empty

    def test_rate_without_significance_not_called(self):
        noisy = {c: ChannelBackground(mean=0.15) for c in CHANNELS}
        sig = _signatures({100: {'mismatch_rate': 0.2}}, coverage=60)
        assert _caller(channel_backgrounds=noisy).call_all('tRNA-x', sig).empty
        # Same site passes on the rate gate alone when testing is disabled
        no_test = _caller(channel_backgrounds=noisy, statistical_test=False)
        assert len(no_test.call_all('tRNA-x', sig)) == 1

    def test_position_one_rt_stop_ignored(self):
        calls = _caller().call_all('tRNA-x', _signatures({1: {'rt_stop_pct': 100.0}}))
        assert calls.empty


class TestBackgroundModel:

    def test_fit_is_calibrated_on_overdispersed_null(self):
        rng = np.random.default_rng(1)
        n_trna, length, cov = 60, 76, 3000
        mu, rho = 1e-3, 1e-2
        a, b = mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho
        pscm, ref = {}, {}
        for i in range(n_trna):
            mat = np.zeros((length, 8))
            stops = rng.binomial(cov, rng.beta(a, b, length))
            mat[:, 2] = cov          # all G
            mat[:, 6] = cov
            mat[:, 7] = stops
            pscm[f't{i}'] = mat
            ref[f't{i}'] = {'seq': 'G' * length}
        bg = estimate_channel_backgrounds(pscm, ref)['rt_stop']
        assert bg.source == 'empirical_bulk'
        assert bg.mean == pytest.approx(mu, rel=0.3)
        assert bg.rho > 1e-3

        k = np.concatenate([m[1:, 7] for m in pscm.values()])
        n = np.full(len(k), cov)
        assert (bg.sf(k, n) < 0.01).mean() < 0.02
        # A single binomial rate is badly anti-conservative on the same data
        assert (ChannelBackground(mean=bg.mean).sf(k, n) < 0.01).mean() > 0.05

    def test_synthetic_tRNAs_preferred(self):
        mat = np.zeros((20, 8))
        mat[:, 2] = mat[:, 6] = 1000
        bgs = estimate_channel_backgrounds(
            {'Synthetic_a': mat, 'native': mat.copy()},
            {'Synthetic_a': {'seq': 'G' * 20}, 'native': {'seq': 'G' * 20}})
        assert all(bg.source == 'synthetic' for bg in bgs.values())

    def test_legacy_single_rate_maps_to_all_channels(self):
        caller = ModificationCaller(background_error_rate=0.02)
        assert {bg.mean for bg in caller.channel_backgrounds.values()} == {0.02}

    def test_larger_fdr_family_is_more_conservative(self):
        p = np.array([1e-4, 1e-3, 1e-2])
        assert (_bh_qvalues(p, 1000) >= _bh_qvalues(p, 3)).all()


class TestPriorsDoNotGateDetection:

    SITES = {100: {'gap_rate': 0.3}, 120: {'rt_stop_pct': 40.0},
             130: {'mismatch_rate': 0.3}}

    def test_same_sites_for_every_enzyme(self):
        sig = _signatures(self.SITES)
        called = {
            enzyme: set(_caller(rt_enzyme=enzyme, rt_temp=55).call_all('t', sig)['position'])
            for enzyme in (None, 'Maxima', 'SSIV', 'TGIRT', 'Indura')
        }
        assert all(sites == {100, 120, 130} for sites in called.values())

    def test_priors_scale_confidence(self):
        # An RT stop is routine for Indura (prior weight 0.74) and suppressed
        # for Maxima at 55 C (0.16), so confidence is scaled accordingly
        sig = _signatures({120: {'rt_stop_pct': 40.0}})
        agnostic = _caller().call_all('t', sig).iloc[0]
        indura = _caller(rt_enzyme='Indura').call_all('t', sig).iloc[0]
        maxima = _caller(rt_enzyme='Maxima', rt_temp=55).call_all('t', sig).iloc[0]
        assert agnostic['channel_prior_factor'] == 1.0
        assert indura['confidence'] > agnostic['confidence'] > maxima['confidence']

    def test_unexpected_channel_flags_signal_the_enzyme_rarely_gives(self):
        # Indura barely deletes (prior weight 0.07); Maxima does (0.35)
        sig = _signatures({100: {'gap_rate': 0.3}})
        indura = _caller(rt_enzyme='Indura').call_all('t', sig).iloc[0]
        maxima = _caller(rt_enzyme='Maxima').call_all('t', sig).iloc[0]
        assert indura['unexpected_channel'] and not maxima['unexpected_channel']
        # flagged, not withheld
        assert indura['modification'] or indura['source']


class TestSiteLevelOutput:

    def test_one_row_per_site_with_ranked_candidates(self):
        calls = _caller().call_all('t', _signatures({46: {'mismatch_rate': 0.3}}))
        assert len(calls) == 1
        assert calls.iloc[0]['n_candidates'] > 1

    def test_sprinzl_typical_position_not_used_as_linear(self):
        # m7G is canonically Sprinzl 46; linear 46 alone must not label it
        row = _caller().call_all('t', _signatures({46: {'mismatch_rate': 0.3}})).iloc[0]
        assert row['modification'] == 'novel_candidate'
        assert not row['in_typical_position']

    def test_unsupported_identity_is_novel_but_still_called(self):
        calls = _caller().call_all('t', _signatures({100: {'gap_rate': 0.3}}))
        row = calls.iloc[0]
        assert row['source'] == 'novel_candidate'
        assert row['identity_support'] == 'ref_nt_only'
        assert row['n_candidates'] > 0

    def test_signature_override_restricts_channel(self):
        # m7G restricted to mismatch: a deletion at 46 can no longer be m7G
        sig = _signatures({46: {'gap_rate': 0.3}})
        default = _caller().call_all('t', sig).iloc[0]
        restricted = _caller(signature_overrides={'m7G': 'mismatch'}).call_all('t', sig).iloc[0]
        assert 'm7G' in default['candidates'].split(';')
        assert 'm7G' not in restricted['candidates'].split(';')

    def test_sample_level_finalize_matches_columns(self):
        caller = _caller()
        sig = _signatures({100: {'gap_rate': 0.3}})
        raw = caller.call_all('t', sig, finalize=False)
        final = caller.finalize_calls(raw, caller.count_tests(sig, 50))
        assert final.equals(caller.call_all('t', sig))
        assert not any(c.startswith(('gate_', 'rate_')) for c in final.columns)


def _known(*pos_names):
    return pd.DataFrame([{'linear_position': p, 'modification_short_name': n}
                         for p, n in pos_names])


class TestSiteIdentity:
    """MODOMICS-first identity at the base the signal implicates."""

    def test_modomics_labels_substitution_site(self):
        row = _caller().call_all('t', _signatures({47: {'mismatch_rate': 0.3}}),
                                 known_mods_df=_known((47, 'acp3U'))).iloc[0]
        assert row['modification'] == 'acp3U'
        assert row['source'] == 'known_modomics'
        assert row['identity_support'] == 'modomics'
        assert row['modified_position'] == 47
        assert row['candidates'].split(';')[0] == 'acp3U'

    def test_rt_stop_implicates_base_5prime_of_stop(self):
        # Stop recorded at 38; m1G is at 37
        row = _caller().call_all('t', _signatures({38: {'rt_stop_pct': 50.0}}),
                                 known_mods_df=_known((37, 'm1G'))).iloc[0]
        assert row['position'] == 38
        assert row['modification'] == 'm1G'
        assert row['modified_position'] == 37

    def test_rt_stop_prefers_preceding_base_when_both_modified(self):
        row = _caller().call_all('t', _signatures({47: {'rt_stop_pct': 50.0}}),
                                 known_mods_df=_known((46, 'm7G'), (47, 'acp3U'))).iloc[0]
        assert row['modification'] == 'm7G'
        assert row['candidates'].split(';')[:2] == ['m7G', 'acp3U']

    def test_deletion_prefers_same_base_then_preceding(self):
        same = _caller().call_all('t', _signatures({47: {'gap_rate': 0.3}}),
                                  known_mods_df=_known((46, 'm7G'), (47, 'acp3U'))).iloc[0]
        prev = _caller().call_all('t', _signatures({47: {'gap_rate': 0.3}}),
                                  known_mods_df=_known((46, 'm7G'))).iloc[0]
        assert same['modification'] == 'acp3U'
        assert prev['modification'] == 'm7G' and prev['modified_position'] == 46

    def test_substitution_does_not_borrow_neighbouring_mod(self):
        row = _caller().call_all('t', _signatures({48: {'mismatch_rate': 0.3}}),
                                 known_mods_df=_known((47, 'acp3U'))).iloc[0]
        assert row['source'] == 'novel_candidate'

    def test_specific_pattern_labels_without_modomics(self):
        # C->T dominant substitution at a C: m3C (C->T, weight 1.2) ranks first
        sig = _signatures({100: {'mismatch_rate': 0.3, 'correct_nt': 'C'}})
        pscm = pd.DataFrame(0, index=range(150), columns=list('ACGTUN-'))
        pscm.loc[99, ['C', 'T']] = [1400, 600]
        row = _caller().call_all('t', sig, pscm_df=pscm).iloc[0]
        assert row['source'] == 'known'
        assert row['identity_support'] == 'pattern'
        assert row['modification'] == 'm3C'

    def test_modomics_overrides_pattern_label(self):
        sig = _signatures({100: {'mismatch_rate': 0.3, 'correct_nt': 'C'}})
        pscm = pd.DataFrame(0, index=range(150), columns=list('ACGTUN-'))
        pscm.loc[99, ['C', 'T']] = [1400, 600]
        row = _caller().call_all('t', sig, pscm_df=pscm,
                                 known_mods_df=_known((100, 'ac4C'))).iloc[0]
        assert row['modification'] == 'ac4C'
        assert 'm3C' in row['candidates'].split(';')

    def test_identity_never_changes_detection(self):
        sig = _signatures({38: {'rt_stop_pct': 50.0}, 100: {'gap_rate': 0.3}})
        with_map = _caller().call_all('t', sig, known_mods_df=_known((37, 'm1G')))
        without = _caller().call_all('t', sig)
        no_priors = _caller(use_position_priors=False).call_all(
            't', sig, known_mods_df=_known((37, 'm1G')))
        assert set(with_map['position']) == set(without['position']) == {38, 100}
        assert (no_priors['source'] == 'novel_candidate').all()


def test_analyzer_mismatch_rate_excludes_deletions():
    analyzer = RTSignatureAnalyzer(min_coverage=10, verbose=False)
    pscm = pd.DataFrame({'A': [0, 0], 'C': [0, 0], 'G': [700, 650], 'T': [0, 50],
                         'U': [0, 0], 'N': [0, 0], '-': [300, 300]})
    mm = analyzer.calculate_mismatch_rates(pscm, 'GG')
    assert mm['mismatch_rate'].tolist() == pytest.approx([0.0, 0.05])
