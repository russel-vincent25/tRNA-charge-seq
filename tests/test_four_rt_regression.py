"""
Regression test on real four-RT data (CHANNEL_FIX_SPEC change 6).

The fixture (tests/data/four_rt_fixture.parquet) holds per-position counts
for 5 E. coli tRNAs x 36 libraries: 4 reverse transcriptases (Indura,
Maxima, SSIV, TGIRT) x 3 temperatures x 3 replicates, from the 2024-06-19
RT comparison (stage 6 mismatch_profile + rt_profile). The tRNAs were
chosen so that between them every enzyme shows signal in every channel:
Arg-ACG (all channels), Lys-TTT (the spec's worked example, position 47),
Phe-GAA (RT-stop dominant), Pro-TGG (deletion dominant), Leu-GAG (no
deletion signal).

The assertions are enzyme-agnostic: they state what must hold whichever
RT produced the data, rather than an ordering of enzymes. "Signal" below
uses the calling thresholds themselves, so the recall test checks that the
whole path -- analyzer, per-sample backgrounds, sample-wide FDR -- keeps
every threshold-passing site; the method is tested by the synthetic
channel-invariance and calibration tests in test_channel_calling.py.

Set TRNASEQ_FOUR_RT_DIR to the full run's modification_analysis directory
to repeat the recall test on all 49 tRNAs.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from trnaseq.modifications.channel_priors import DERIVATION_THRESHOLDS as THRESHOLDS
from trnaseq.modifications.modification_caller import (
    CHANNELS,
    ChannelBackground,
    ModificationCaller,
    estimate_channel_backgrounds,
)
from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer

FIXTURE = Path(__file__).parent / 'data' / 'four_rt_fixture.parquet'
MIN_COVERAGE = 100
COUNT_COLS = ['A_count', 'C_count', 'G_count', 'T_count', 'N_count',
              'gap_count', 'coverage', 'rt_stop_count']


def _channel_rates(df):
    """Per-position substitution / deletion / RT-stop rates from raw counts."""
    correct = np.choose(df['ref_nt'].map({'A': 0, 'C': 1, 'G': 2, 'T': 3}).to_numpy(),
                        df[['A_count', 'C_count', 'G_count', 'T_count']].to_numpy().T)
    cov = df['coverage'].to_numpy(float)
    return pd.DataFrame({
        'mismatch': (cov - correct - df['gap_count']) / cov,
        'deletion': df['gap_count'] / cov,
        'rt_stop': np.where(df['position'] > 1, df['rt_stop_count'] / cov, 0.0),
    }, index=df.index)


def _call_sample(df, rt_enzyme=None):
    """Stage 6 calling path for one sample's rows."""
    pscm, ref = {}, {}
    for name, g in df.sort_values(['tRNA_name', 'position']).groupby('tRNA_name', observed=True):
        pscm[name] = g[COUNT_COLS].to_numpy(float)
        ref[name] = {'seq': ''.join(g['ref_nt'].astype(str)), 'seq_len': len(g)}
    analyzer = RTSignatureAnalyzer(min_coverage=MIN_COVERAGE, verbose=False)
    analyzer.reference_sequences = ref
    caller = ModificationCaller(
        organism='Escherichia coli', rt_enzyme=rt_enzyme,
        channel_backgrounds=estimate_channel_backgrounds(pscm, ref, min_coverage=MIN_COVERAGE))
    calls, n_tests = [], 0
    for name, pscm_df in analyzer.load_pscm_from_positional(pscm).items():
        sig = analyzer.analyze_trna_with_actual_stops(
            name, pscm_df, rt_stop_counts=pscm[name][:, 7])['signatures']
        c = caller.call_all(name, sig, pscm_df, ref[name]['seq'],
                            min_coverage=MIN_COVERAGE, finalize=False)
        n_tests += caller.count_tests(sig, MIN_COVERAGE)
        if not c.empty:
            calls.append(c)
    if not calls:
        return pd.DataFrame(columns=['trna_name', 'position', 'channels_fired'])
    return caller.finalize_calls(pd.concat(calls, ignore_index=True), n_tests)


def _call_all_samples(fixture):
    out = []
    for sample, df in fixture.groupby('sample', observed=True):
        calls = _call_sample(df, rt_enzyme=str(df['rt_enzyme'].iloc[0]))
        out.append(calls.assign(sample=sample, rt_enzyme=str(df['rt_enzyme'].iloc[0])))
    return pd.concat(out, ignore_index=True)


def _signal_sites(fixture):
    """Site-observations at or above any channel threshold (coverage >= 100)."""
    f = fixture[fixture['coverage'] >= MIN_COVERAGE]
    rates = _channel_rates(f)
    hit = pd.DataFrame({c: rates[c] >= THRESHOLDS[c] for c in CHANNELS})
    sig = f.loc[hit.any(axis=1), ['sample', 'rt_enzyme', 'tRNA_name', 'position']].copy()
    sig['mismatch_only'] = hit.loc[sig.index, 'mismatch'].to_numpy()
    sig[['sample', 'rt_enzyme', 'tRNA_name']] = sig[['sample', 'rt_enzyme', 'tRNA_name']].astype(str)
    return sig


def _recall(signal, calls):
    called = set(zip(calls['sample'], calls['trna_name'], calls['position']))
    signal = signal.assign(called=[k in called for k in zip(
        signal['sample'], signal['tRNA_name'], signal['position'])])
    return signal.groupby('rt_enzyme').agg(
        new=('called', 'mean'), mismatch_only=('mismatch_only', 'mean'), n=('called', 'size'))


@pytest.fixture(scope='module')
def fixture():
    return pd.read_parquet(FIXTURE)


@pytest.fixture(scope='module')
def calls(fixture):
    return _call_all_samples(fixture)


def test_fixture_exercises_the_problem(fixture):
    """Mismatch-only recall must be uneven across enzymes, or the
    recall test below would pass without testing anything."""
    r = _recall(_signal_sites(fixture), pd.DataFrame(columns=['sample', 'trna_name', 'position']))
    assert len(r) == 4 and (r['n'] >= 50).all()
    assert r['mismatch_only'].max() - r['mismatch_only'].min() > 0.15
    assert r['mismatch_only'].max() < 0.6


def test_recall_even_across_enzymes(fixture, calls):
    """Every enzyme's signal sites are called, whichever channel carries
    them, so recall no longer depends on the RT."""
    r = _recall(_signal_sites(fixture), calls)
    assert (r['new'] >= 0.95).all(), r
    assert r['new'].max() - r['new'].min() <= 0.05, r
    assert (r['new'] > r['mismatch_only']).all(), r


def test_no_regression_on_substitution_sites(fixture, calls):
    """Every site a mismatch-only caller would see is still called, and
    still attributed to the substitution channel."""
    sig = _signal_sites(fixture)
    sub = sig[sig['mismatch_only']]
    called = calls.set_index(['sample', 'trna_name', 'position'])['channels_fired']
    keys = list(zip(sub['sample'], sub['tRNA_name'], sub['position']))
    missing = [k for k in keys if k not in called.index]
    assert not missing, missing[:5]
    assert all('mismatch' in called.loc[k] for k in keys)


def test_backgrounds_calibrated_on_real_null(fixture):
    """On each library's null bulk, the fitted per-channel background gives
    about alpha false positives; a single binomial rate gives far more."""
    fitted, binomial = [], []
    for _, df in fixture.groupby('sample', observed=True):
        pscm, ref = {}, {}
        for name, g in df.sort_values('position').groupby('tRNA_name', observed=True):
            pscm[name] = g[COUNT_COLS].to_numpy(float)
            ref[name] = {'seq': ''.join(g['ref_nt'].astype(str))}
        bgs = estimate_channel_backgrounds(pscm, ref, min_coverage=MIN_COVERAGE)
        f = df[(df['coverage'] >= MIN_COVERAGE) & (df['position'] > 1)]
        rates = _channel_rates(f)
        n = f['coverage'].to_numpy()
        for c in ('mismatch', 'rt_stop'):
            bulk = (rates[c] < 0.5 * THRESHOLDS[c]).to_numpy()
            k = np.rint(rates[c].to_numpy() * n)[bulk]
            fitted.append((bgs[c].sf(k, n[bulk]) < 0.01).mean())
            binomial.append((ChannelBackground(mean=bgs[c].mean).sf(k, n[bulk]) < 0.01).mean())
    assert np.median(fitted) <= 0.02
    assert np.median(binomial) > 2 * np.median(fitted)


def test_enzyme_label_never_changes_detection(fixture):
    """Channel priors scale confidence only: the same library called with
    and without its RT recorded yields the same sites."""
    for sample in ('Indura-60-R1', 'Maxima-55-R1'):
        df = fixture[fixture['sample'] == sample]
        with_rt = _call_sample(df, rt_enzyme=str(df['rt_enzyme'].iloc[0]))
        without = _call_sample(df, rt_enzyme=None)
        assert set(zip(with_rt['trna_name'], with_rt['position'])) == \
            set(zip(without['trna_name'], without['position']))


@pytest.mark.skipif(not os.environ.get('TRNASEQ_FOUR_RT_DIR'),
                    reason='set TRNASEQ_FOUR_RT_DIR to run on the full four-RT dataset')
def test_recall_even_across_enzymes_full_dataset():
    rows = []
    for d in sorted(Path(os.environ['TRNASEQ_FOUR_RT_DIR']).glob('*-*-R*')):
        mm = pd.read_parquet(d / 'mismatch_profile.parquet')
        rp = pd.read_parquet(d / 'rt_profile.parquet')
        m = mm.merge(rp[['tRNA_name', 'position', 'rt_stop_count']], on=['tRNA_name', 'position'])
        rows.append(m.assign(sample=d.name, rt_enzyme=d.name.split('-')[0]))
    full = pd.concat(rows, ignore_index=True)
    r = _recall(_signal_sites(full), _call_all_samples(full))
    assert (r['new'] >= 0.95).all(), r
    assert r['new'].max() - r['new'].min() <= 0.05, r
