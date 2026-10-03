"""
Tests for RT-enzyme channel handling (CHANNEL_FIX_SPEC changes 1-3):
- rt_enzyme / rt_temp provenance on the caller and its calls
- 'combined' default signature_type with per-modification override
- Per-enzyme channel priors: derivation and lookup
"""

import warnings

import numpy as np
import pandas as pd
import pytest

from trnaseq.modifications.channel_priors import (
    CHANNELS,
    derive_channel_priors,
    load_channel_priors,
    normalize_rt_enzyme,
)
from trnaseq.modifications.modification_caller import (
    MODIFICATION_PROFILES,
    ModificationCaller,
)


# ---------------------------------------------------------------------------
# Change 1: provenance
# ---------------------------------------------------------------------------

class TestRTProvenance:

    def test_defaults_are_null_safe(self):
        caller = ModificationCaller()
        assert caller.rt_enzyme is None
        assert caller.rt_temp is None
        assert caller.channel_priors is None

    def test_calls_record_enzyme_and_temp(self):
        sig = pd.DataFrame([{
            'position': 1, 'has_signature': True, 'mismatch_rate': 0.4,
            'rt_stop_pct': 0.0, 'gap_rate': 0.0, 'coverage': 1000,
        }])
        caller = ModificationCaller(min_confidence=0.0, rt_enzyme='Maxima',
                                    rt_temp=55)
        calls = caller.call_all('tRNA-x', sig)
        assert not calls.empty
        assert (calls['rt_enzyme'] == 'Maxima').all()
        assert (calls['rt_temp'] == 55).all()


# ---------------------------------------------------------------------------
# Change 2: 'combined' default, overridable
# ---------------------------------------------------------------------------

class TestSignatureTypeDefault:

    def test_every_profile_defaults_to_combined(self):
        assert {p.signature_type for p in MODIFICATION_PROFILES.values()} == {'combined'}

    def test_override_applies_without_mutating_defaults(self):
        caller = ModificationCaller(signature_overrides={'m3C': 'mismatch'})
        assert caller.profiles['m3C'].signature_type == 'mismatch'
        assert MODIFICATION_PROFILES['m3C'].signature_type == 'combined'
        assert ModificationCaller().profiles['m3C'].signature_type == 'combined'

    @pytest.mark.parametrize('overrides', [
        {'not_a_mod': 'mismatch'},
        {'m3C': 'deletion'},
    ])
    def test_invalid_override_raises(self, overrides):
        with pytest.raises(ValueError):
            ModificationCaller(signature_overrides=overrides)


# ---------------------------------------------------------------------------
# Change 3: per-enzyme priors
# ---------------------------------------------------------------------------

class TestChannelPriors:

    @pytest.mark.parametrize('name,expected', [
        ('Maxima H Minus', 'Maxima'),
        ('SuperScript IV', 'SSIV'),
        ('TGIRT-III', 'TGIRT'),
        ('Induro', 'Indura'),
        ('AMV', None),
        (None, None),
    ])
    def test_normalize_rt_enzyme(self, name, expected):
        assert normalize_rt_enzyme(name) == expected

    def test_unknown_enzyme_falls_back(self):
        assert load_channel_priors(None) is None
        with pytest.warns(UserWarning, match='Falling back'):
            assert load_channel_priors('AMV') is None

    def test_shipped_priors_cover_four_enzymes(self):
        for enzyme in ('Maxima', 'SSIV', 'TGIRT', 'Indura'):
            priors = load_channel_priors(enzyme)
            assert priors['rt_temp'] is None  # pooled row
            assert sum(priors['weights'].values()) == pytest.approx(1.0, abs=1e-3)

    def test_shipped_priors_dominant_channels(self):
        """Pooled priors reflect the four-RT comparison: Indura RT-stop
        dominant, Maxima deletion-weighted, TGIRT mismatch-weighted."""
        def top(enzyme):
            w = load_channel_priors(enzyme)['weights']
            return max(w, key=w.get)
        assert top('Indura') == 'rt_stop'
        assert top('Maxima') == 'deletion'
        assert top('TGIRT') == 'mismatch'

    def test_temperature_row_and_pooled_fallback(self):
        assert load_channel_priors('Maxima', 55)['rt_temp'] == 55.0
        assert load_channel_priors('Maxima', 42)['rt_temp'] is None

    def test_derive_from_observations(self):
        # Two enzymes, one temp; X sees signal only as deletion, Y only as
        # mismatch; background rows carry no signal and must be ignored.
        obs = pd.DataFrame({
            'rt_enzyme': ['X'] * 3 + ['Y'] * 3,
            'rt_temp': [55.0] * 6,
            'mismatch': [0.01, 0.01, 0.00, 0.30, 0.25, 0.00],
            'deletion': [0.40, 0.20, 0.00, 0.00, 0.00, 0.00],
            'rt_stop': [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
        })
        priors = derive_channel_priors(obs)
        pooled = priors[priors['rt_temp'].isna()].set_index('rt_enzyme')
        assert pooled.loc['X', 'n_signal_obs'] == 2
        assert pooled.loc['X', 'weight_deletion'] == pytest.approx(1.0)
        assert pooled.loc['Y', 'weight_mismatch'] == pytest.approx(1.0)
        assert len(priors) == 4  # (X, Y) x (55, pooled)

    def test_caller_attaches_priors(self):
        caller = ModificationCaller(rt_enzyme='TGIRT-III', rt_temp=60)
        assert caller.channel_priors['rt_enzyme'] == 'TGIRT'
        assert caller.channel_priors['rt_temp'] == 60.0
        assert set(caller.channel_priors['weights']) == set(CHANNELS)
