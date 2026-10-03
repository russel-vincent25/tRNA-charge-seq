"""
Tests for MODOMICS lookup: reference-name keys, the symbol table, and
isotype-donor anticodon-loop handling.
"""

import pandas as pd
import pytest

from trnaseq.modifications.modification_caller import ChannelBackground, ModificationCaller
from trnaseq.modifications.modomics import (
    MODOMICS_SYMBOL_NAMES,
    MODOMICS_TO_BASE,
    MODOMICSAnnotator,
    _lookup_key,
    _strip_modomics_sequence,
)

ECOLI = 'Escherichia_coli_str_K12_substr_MG1655_tRNA-'


@pytest.mark.parametrize('name,key', [
    (ECOLI + 'Ile2-CAT-1-1', ('ile', 'CAT')),
    (ECOLI + 'fMet-CAT-1-1', ('ini', 'CAT')),
    ('Homo_sapiens_tRNA-iMet-CAT-1-1', ('ini', 'CAT')),
    ('Homo_sapiens_tRNA-Ala-AGC-1-1', ('ala', 'AGC')),
    ('Homo_sapiens_mito_tRNA-Leu1-TAG', None),
    ('no_hyphens', None),
])
def test_lookup_key(name, key):
    assert _lookup_key(name) == key


@pytest.mark.parametrize('symbol,name,parent', [
    ('}', 'k2C', 'C'),     # lysidine: a modified C, previously mapped to U
    ('?', 'm5C', 'C'),     # previously 'xG' on G
    ('O', 'm1I', 'A'),
    ('Ч', 'i6A', 'A'),
    ('P', 'Psi', 'U'),     # curated spelling kept
    ('#', 'Gm', 'G'),      # curated; absent from the MODOMICS table
])
def test_symbol_table(symbol, name, parent):
    assert MODOMICS_SYMBOL_NAMES[symbol] == name
    assert MODOMICS_TO_BASE[symbol] == parent


def test_unknown_nucleotide_is_not_a_modification():
    base, mods = _strip_modomics_sequence('AC.G')
    assert base == 'ACNG'
    assert mods == []


@pytest.fixture(scope='module')
def ecoli():
    ann = MODOMICSAnnotator('Escherichia coli')
    ann.get_modifications(use_api=False)
    return ann


def _ref(ann, key):
    base, _ = _strip_modomics_sequence(ann._modomics_sequences[key])
    return base.replace('U', 'T')


def test_ile2_maps_to_its_modomics_sequence(ecoli):
    known = ecoli.get_known_mods_linear(ECOLI + 'Ile2-CAT-1-1', _ref(ecoli, ('ile', 'CAT')))
    assert set(known['mapping']) == {'exact'}
    assert 'k2C' in set(known['modification_short_name'])


def test_mito_reference_gets_no_map(ecoli):
    assert ecoli.get_known_mods_linear(
        'Homo_sapiens_mito_tRNA-Ile-CAT', _ref(ecoli, ('ile', 'CAT'))).empty


def test_donor_anticodon_loop_opt_in(ecoli):
    # Pro-GGG has no MODOMICS sequence of its own; Pro-CGG is the donor
    ref = _ref(ecoli, ('pro', 'CGG'))
    name = ECOLI + 'Pro-GGG-1-1'
    default = ecoli.get_known_mods_linear(name, ref)
    kept = ecoli.get_known_mods_linear(name, ref, include_donor_anticodon_loop=True)
    assert set(default['mapping']) == {'isotype'}
    assert 'm1G' not in set(default['modification_short_name'])
    assert 'm1G' in set(kept['modification_short_name'])


def test_caller_marks_donor_evidence_weaker():
    bgs = {c: ChannelBackground(mean=1e-3) for c in ('mismatch', 'deletion', 'rt_stop')}
    sig = pd.DataFrame([{'position': 38, 'coverage': 2000, 'correct_nt': 'G',
                         'mismatch_rate': 0.3, 'gap_rate': 0.0, 'rt_stop_pct': 0.0}])
    caller = ModificationCaller(channel_backgrounds=bgs)

    def call(mapping):
        known = pd.DataFrame([{'linear_position': 38, 'modification_short_name': 'm1G',
                               'mapping': mapping}])
        return caller.call_all('t', sig, known_mods_df=known).iloc[0]

    exact, donor = call('exact'), call('isotype')
    assert exact['identity_support'] == 'modomics'
    assert donor['identity_support'] == 'modomics_isotype'
    assert donor['modification'] == 'm1G'
    assert exact['confidence'] > donor['confidence']
