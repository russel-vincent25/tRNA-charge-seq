#!/usr/bin/env python3
"""PROTOTYPE — predict the identity of novel modification sites.

Not part of the pipeline. Trains a classifier on sites the caller labelled
from MODOMICS and predicts labels for ``novel_candidate`` sites, using each
site's RT fingerprint: median substitution / deletion / RT-stop rate per
enzyme x temperature condition, the substitution spectrum, and the local
sequence. Validation holds out whole isotypes so the model cannot memorise
sites.

Result on the 2024-06-19 four-RT run (E. coli, 100 labelled sites, 5
classes: D, m7G, ms2i6A, acp3U, m1G), held-out by isotype:

    majority-class baseline                     0.39
    RT rates only, one condition (3 cols)       0.62
    RT rates only, all 12 conditions (36)       0.75
    RT rates + substitution spectrum            0.83
    position only (rel_ac, rel_3p)              0.71
    reference base only                         0.75
    position + reference base                   0.87
    sequence + position, NO rate features       0.88
    everything                                  0.88

Two conclusions, and the second is why this is not in the pipeline.

1. The RT signal does carry identity information -- 0.75 from 36 rate
   numbers against a 0.39 baseline -- and the multi-condition fingerprint
   is worth +0.13 over a single condition (0.75 vs 0.62). That is a
   measured argument for running several RTs rather than one.

2. For these five classes the signal is nonetheless redundant: sequence
   and position alone reach 0.88, and adding every rate feature moves
   nothing. These modifications sit at stereotyped positions (D in the
   D-loop, m1G37, m7G46, acp3U47, ms2i6A37), so a model that knows where
   it is can name them without looking at the RT data.

That makes the predictor close to useless for its stated purpose. At a
genuinely novel site -- no annotation -- the positional features would
predict whatever modification usually sits there, which is circular, and
the part that could generalise (the rates) is the weaker part. An earlier
version of this docstring reported "0.84 from the RT fingerprint"; that
figure mixed rate and sequence features and overstated the signal's
contribution.

Open questions before this goes into the pipeline (see
docs/CHANNEL_COMBINED_DESIGN.md, "Novel-site prediction"):
- accuracy from a single enzyme/condition (normal runs) rather than the
  12-condition fingerprint
- other organisms and reference sets
- ground truth beyond MODOMICS (e.g. modification-enzyme knockouts)

Usage:
    python docs/prototypes/novel_site_classifier.py \\
        --profiles <run>/modification_analysis \\
        --organism "Escherichia coli" -o novel_site_predictions.csv

``--profiles`` holds one directory per sample (named <enzyme>-<temp>-<rep>)
with mismatch_profile.parquet and rt_profile.parquet from stage 6.
"""

import argparse
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold

from trnaseq.modifications.modification_caller import (
    ModificationCaller,
    estimate_channel_backgrounds,
)
from trnaseq.modifications.modomics import MODOMICSAnnotator
from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer

warnings.filterwarnings('ignore')
PSCM_COLS = ['A_count', 'C_count', 'G_count', 'T_count', 'N_count',
             'gap_count', 'coverage', 'rt_stop_count']


def load_profiles(directory, pattern):
    """{sample: merged per-position profile} plus (enzyme, temp) per sample."""
    regex = re.compile(pattern)
    out, meta = {}, {}
    for d in sorted(Path(directory).iterdir()):
        m = regex.match(d.name)
        if not m or not (d / 'mismatch_profile.parquet').exists():
            continue
        mm = pd.read_parquet(d / 'mismatch_profile.parquet')
        rp = pd.read_parquet(d / 'rt_profile.parquet')
        out[d.name] = mm.merge(rp[['tRNA_name', 'position', 'rt_stop_count',
                                   'rt_stop_fraction']], on=['tRNA_name', 'position'])
        meta[d.name] = (m.group('rt_enzyme'), m.group('rt_temp'))
    return out, meta


def call_sample(df, annotator, known_cache, rt_enzyme, rt_temp):
    """Re-call one sample with the current caller (labels come from here)."""
    pscm, ref = {}, {}
    for t, g in df.sort_values(['tRNA_name', 'position']).groupby('tRNA_name'):
        pscm[t] = g[PSCM_COLS].to_numpy(float)
        ref[t] = {'seq': ''.join(g['ref_nt']), 'seq_len': len(g)}
    analyzer = RTSignatureAnalyzer(min_coverage=100, verbose=False)
    analyzer.reference_sequences = ref
    dfs = analyzer.load_pscm_from_positional(pscm)
    caller = ModificationCaller(
        rt_enzyme=rt_enzyme, rt_temp=float(rt_temp),
        channel_backgrounds=estimate_channel_backgrounds(pscm, ref, min_coverage=100))
    calls, n_tests = [], 0
    for t, pdf in dfs.items():
        if t not in known_cache:
            known_cache[t] = annotator.get_known_mods_linear(
                t, ref[t]['seq'], include_donor_anticodon_loop=True)
        sig = analyzer.analyze_trna_with_actual_stops(t, pdf, rt_stop_counts=pscm[t][:, 7])['signatures']
        c = caller.call_all(t, sig, pdf, ref[t]['seq'], min_coverage=100,
                            known_mods_df=known_cache[t], finalize=False)
        n_tests += caller.count_tests(sig, 100)
        if not c.empty:
            calls.append(c)
    return caller.finalize_calls(pd.concat(calls, ignore_index=True), n_tests) if calls else pd.DataFrame()


def anticodon_start(name, seq):
    """1-based start of the anticodon: the match to the name's anticodon nearest
    linear 34 within 28-42 (tRNA lengths vary, so no fixed offset)."""
    ac = name.split('-')[2].upper().replace('U', 'T')
    hits = [i + 1 for i in range(27, 42) if seq[i:i + 3] == ac]
    return min(hits, key=lambda h: abs(h - 34)) if hits else np.nan


def build_features(profiles, meta):
    rows = []
    for sample, df in profiles.items():
        enzyme, temp = meta[sample]
        rows.append(df.assign(cond=f'{enzyme}|{temp}',
                              stop=np.where(df.position > 1, df.rt_stop_fraction, 0.0)))
    raw = pd.concat(rows, ignore_index=True)
    ref = {t: ''.join(g.drop_duplicates('position').sort_values('position').ref_nt)
           for t, g in raw.groupby('tRNA_name')}

    fp = raw.pivot_table(index=['tRNA_name', 'position'], columns='cond',
                         values=['mismatch_rate', 'deletion_rate', 'stop'], aggfunc='median')
    fp.columns = [f'{a}|{b}' for a, b in fp.columns]
    spec = raw.groupby(['tRNA_name', 'position'])[['A_count', 'C_count', 'G_count', 'T_count']].sum()
    X = fp.join(spec).reset_index()

    X['ref'] = [ref[t][p - 1] for t, p in zip(X.tRNA_name, X.position)]
    X['prev'] = [ref[t][p - 2] if p > 1 else 'N' for t, p in zip(X.tRNA_name, X.position)]
    X['next'] = [ref[t][p] if p < len(ref[t]) else 'N' for t, p in zip(X.tRNA_name, X.position)]
    counts = X[['A_count', 'C_count', 'G_count', 'T_count']].to_numpy(float)
    is_ref = np.column_stack([(X.ref == nt).to_numpy() for nt in 'ACGT'])
    nonref = np.where(is_ref, 0, counts)
    for i, nt in enumerate('ACGT'):
        X[f'spec_{nt}'] = nonref[:, i] / np.maximum(nonref.sum(1), 1)
    X['rel_ac'] = X.position - np.array([anticodon_start(t, ref[t]) for t in X.tRNA_name])
    X['rel_3p'] = X.position - np.array([len(ref[t]) for t in X.tRNA_name])
    for col in ('ref', 'prev', 'next'):
        X = X.join(pd.get_dummies(X[col], prefix=col).astype(float))

    signature = [c for c in X.columns if '|' in c or c.startswith(('spec_', 'ref_', 'prev_', 'next_'))]
    X[signature + ['rel_ac', 'rel_3p']] = X[signature + ['rel_ac', 'rel_3p']].fillna(0)
    return X, signature


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--profiles', required=True)
    ap.add_argument('--organism', default='Escherichia coli')
    ap.add_argument('--sample-pattern', default=r'^(?P<rt_enzyme>[^-]+)-(?P<rt_temp>\d+)-')
    ap.add_argument('--min-class-sites', type=int, default=4)
    ap.add_argument('-o', '--output', default='novel_site_predictions.csv')
    args = ap.parse_args(argv)

    profiles, meta = load_profiles(args.profiles, args.sample_pattern)
    print(f'{len(profiles)} samples, conditions: {sorted(set(meta.values()))}')

    annotator = MODOMICSAnnotator(args.organism)
    annotator.get_modifications(use_api=False)
    known_cache = {}
    calls = pd.concat([call_sample(df, annotator, known_cache, *meta[s]).assign(sample=s)
                       for s, df in profiles.items()], ignore_index=True)

    # Labels: exact-isodecoder MODOMICS calls only (borrowed maps are weaker)
    exact = calls[calls.identity_support == 'modomics']
    labels = exact.groupby(['trna_name', 'position']).modification.agg(lambda s: s.mode()[0])
    novel = calls[calls.source == 'novel_candidate'].groupby(['trna_name', 'position']).size()

    X, signature = build_features(profiles, meta)
    key = pd.MultiIndex.from_frame(X[['tRNA_name', 'position']])
    X['label'] = labels.reindex(key).values
    X['novel_obs'] = novel.reindex(key).values

    L = X[X.label.notna()]
    counts = L.label.value_counts()
    L = L[L.label.isin(counts[counts >= args.min_class_sites].index)]
    print(f'labelled sites: {len(L)}; classes: {L.label.value_counts().to_dict()}')
    groups = L.tRNA_name.str.extract(r'tRNA-([A-Za-z]+)')[0]
    n_splits = min(5, groups.nunique())

    def held_out_accuracy(cols):
        pred = pd.Series(index=L.index, dtype=object)
        for tr, te in GroupKFold(n_splits=n_splits).split(L, L.label, groups):
            m = RandomForestClassifier(n_estimators=500, class_weight='balanced', random_state=0)
            pred.iloc[te] = m.fit(L.iloc[tr][cols], L.iloc[tr].label).predict(L.iloc[te][cols])
        return (pred == L.label).mean()

    print(f"{'majority class':<36}{L.label.value_counts(normalize=True).iloc[0]:.2f}")
    print(f"{'reference base only':<36}{held_out_accuracy([c for c in signature if c.startswith('ref_')]):.2f}")
    print(f"{'RT fingerprint + sequence':<36}{held_out_accuracy(signature):.2f}")
    full = signature + ['rel_ac', 'rel_3p']
    print(f"{'+ anticodon-relative position':<36}{held_out_accuracy(full):.2f}")

    model = RandomForestClassifier(n_estimators=500, class_weight='balanced', random_state=0)
    model.fit(L[full], L.label)
    N = X[X.novel_obs.notna()].copy()
    if N.empty:
        print('no novel sites to predict')
        return
    proba = model.predict_proba(N[full])
    N['predicted'] = model.classes_[proba.argmax(1)]
    N['probability'] = proba.max(1).round(3)
    out = N[['tRNA_name', 'position', 'ref', 'novel_obs', 'predicted', 'probability']]
    out.sort_values('novel_obs', ascending=False).to_csv(args.output, index=False)
    print(f'{len(out)} novel sites -> {args.output}')


if __name__ == '__main__':
    main()
