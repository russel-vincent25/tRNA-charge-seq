# Design: channel-aware calling and site identity (CHANNEL_FIX_SPEC change 4)

**Status (2026-10-03):** detection (change 4, §1–§5), site identity (change 4b, §6–§7) and the
novel-site lookup fixes and table (§8) are implemented. The novel-site predictor (§9) is a
prototype. The N-masked-reference fix (§8), change 5 (§10) and change 6 (§11) are done; every item in the spec's "Done when" list is met.

Evidence comes from the 2024-06-19 four-RT run (4 enzymes × 3 temperatures × 3 replicates,
`ecoli.fa`, 49 references). Only `mismatch_profile.parquet` and `rt_profile.parquet` were used:
110,114 site-observations at coverage ≥ 100, with position 1 excluded.

Decisions taken:

| | Question | Decision |
|---|---|---|
| A | How channels combine | **OR** across substitution, deletion and RT stop |
| B | Role of per-enzyme channel priors | **Confidence and QC flag only**, never detection |
| C | Output shape | **One row per site** with a ranked candidate list |
| E | Fix site identity before change 5 | **Yes** (change 4b, separate commit) |
| F | Source of position evidence | **MODOMICS**, mapped by sequence alignment; no Sprinzl table. MODOMICS carries no Sprinzl numbering (the API returns only `seq`, `subtype`, `anticodon`, `organism`, `type` and empty ID fields) |
| G | Commit layout | 4b as its own commit, so it can be reverted alone |
| D | Regression test (change 6) | **Enzyme-agnostic**: channel invariance, calibration, no-regression, even recall across enzymes. Not "Indura gains most / TGIRT least". |

---

## 1. What `combined` did before this change

It did not OR the channels, require consensus, or sum evidence. **It was the mismatch path with a
confidence bonus.**

- `call_modification_at_position` gated every profile on `mismatch_rate >= profile.min_rate`,
  whatever its `signature_type`. `'rt_stop'` and `'gap'` were accepted values with no code path.
- `combined` only added an RT-stop term to `calculate_confidence`. So change 2 (default everything to
  `combined`) recovered no signal on its own.
- **The caller's "mismatch" included deletions.** `RTSignatureAnalyzer.calculate_mismatch_rates` used
  `(coverage − correct)/coverage`, with gaps in the coverage. A deletion-dominated site therefore
  passed the rate gate. It then failed any specific substitution pattern, because the gaps sat in the
  denominator, and fell through to `novel_candidate`.
- **0.001 was never a hardcoded constant.** It is the `FLOOR` in `estimate_background_error_rate`.
  Without synthetic spike-ins, that function takes the 25th percentile of per-position rates. That
  percentile is 0, because 41% of positions have zero mismatches, so the floor wins.
- **`min_confidence = 0.5` was a hidden detection gate.** `tests/benchmarks/modification_signal_audit.py`
  shows it implies a mismatch-rate floor of about 0.31 at typical positions.

## 2. Background model: per channel, per sample, overdispersed

A single binomial rate is not calibrated for any channel. A beta-binomial fitted by maximum
likelihood, per sample and per channel, is. The fit uses null-bulk positions only: rate below half the
channel threshold, at most 5,000 positions, seeded subsample. Values are medians over each enzyme's 9
samples:

| RT | null mean mm / del / stop | overdispersion ρ mm / del / stop | null positions at p<0.01: binomial → beta-binomial |
|---|---|---|---|
| Indura | 8.2e-4 / 0.7e-4 / 10.6e-4 | 0.0011 / 0.0027 / 0.0141 | mm 7.0%→1.7%, stop 4.5%→1.2% |
| Maxima | 8.5e-4 / 2.1e-4 / 8.3e-4 | 0.0012 / 0.0058 / 0.0103 | mm 5.3%→1.0%, stop 4.5%→1.2% |
| SSIV | 8.2e-4 / 1.6e-4 / 9.9e-4 | 0.0012 / 0.0050 / 0.0116 | mm 6.5%→1.0%, stop 4.2%→1.2% |
| TGIRT | 7.6e-4 / 1.2e-4 / 9.8e-4 | 0.0008 / 0.0053 / 0.0111 | mm 6.7%→1.2%, stop 5.3%→1.3% |

The RT-stop null is about 10× more overdispersed than mismatch. A per-channel *rate* plugged into a
binomial test would stay miscalibrated, worst in exactly the channel Indura depends on. Synthetic
spike-ins, when present, feed the same fit.

## 3. Aggregation rules compared

A channel *fires* when its rate is at or above the channel threshold (mismatch 0.10, deletion 0.10,
RT stop 0.20; unchanged) **and** its BH q < 0.01. The false-positive proxy is replicate
reproducibility: of the site×temperature units called in at least one replicate, the share called
in all three.

| rule | site-obs called (Indura / Maxima / SSIV / TGIRT) | reproduced 3/3 |
|---|---|---|
| mismatch only (old behaviour) | 82 / 130 / 84 / 137 | 0.86 / 0.73 / 0.69 / 0.68 |
| **OR** | 282 / 268 / 237 / 272 | 0.67 / 0.62 / 0.53 / 0.52 |
| ≥ 2 channels | 20 / 71 / 52 / 66 | 0.86 / 0.71 / 0.50 / 0.39 |
| Fisher sum | 285 / 280 / 260 / 284 | 0.66 / 0.62 / 0.50 / 0.52 |

- **Requiring ≥ 2 channels is ruled out.** Indura's signal is single-channel, so it would lose 93%.
- **OR is identical to OR-ing the thresholds at this depth.** The thresholds decide calls; the null
  model decides p-values.
- **A Fisher sum adds 1–10%** at the same reproducibility. It also assumes independence, which
  substitution and deletion don't have (they compete for the same reads).
- **OR's extra calls are less reproducible.** This is the measurable cost of OR. Part of it is
  near-threshold depth noise, which I can't separate from false positives without ground truth.
- **Gains over mismatch-only:** Indura ×3.4, SSIV ×2.8, Maxima ×2.06, TGIRT ×1.99. The 3% gap
  between Maxima and TGIRT is why decision D drops the ordering assertion.

## 4. Implemented (change 4)

`trnaseq/modifications/modification_caller.py`:

- **`ChannelBackground` and `estimate_channel_backgrounds()`:** a beta-binomial null per channel, fitted
  per sample. Synthetic tRNAs are preferred; otherwise the empirical bulk is used. RT stops at
  position 1 are excluded (every full-length read ends there).
- **Detection:** OR over channels. Each channel needs rate ≥ threshold, plus BH q < α over all
  site × channel tests **in the sample**. The pipeline gets this through `call_all(finalize=False)`,
  then `count_tests()`, then `finalize_calls()`. Tests not reported count as larger p-values, which is
  conservative.
- **Substitution-only mismatch:** `calculate_mismatch_rates` subtracts gaps, and
  `match_mismatch_pattern` excludes `-` and `N` from the substitution total.
- **One row per site**, with these columns:
  - per-channel `pvalue_*`, `qvalue_*`, `fold_change_*`
  - `channels_fired`, `n_channels_fired`, `dominant_channel` (largest rate/threshold among fired
    channels)
  - `pvalue`: Bonferroni over channels, for Fisher combination across replicates
  - `fold_change` and `background_error_rate`: those of the dominant channel
  - `candidates` (ranked), `n_candidates`, `identity_support`
- **Priors (decision B):** `confidence` is multiplied by `1 + (w_fired_max − 1/3)`. `unexpected_channel`
  is set when every fired channel has prior weight < 0.10. With no enzyme recorded the factor is 1.0.
  The same 3,115 sites were called on the example project with TGIRT priors and with none.
- **`signature_type`** now restricts which channels can support a profile as a candidate.
- **`min_confidence` defaults to 0.** Confidence ranks calls; it does not gate them.
- **`discover_novel` is deprecated and ignored.** Sites without supported identity are always reported.
  Filtering them would drop most detected sites under the default config.
- **Pipeline and CLI:** stage 6 Phase 2 fits backgrounds per sample and logs per-channel medians.
  `modification_summary.csv` gains `bg_mean_*`, `bg_rho_*`, `rt_enzyme` and `rt_temp`. The standalone
  CLI does the same. Previously it ran a fixed 0.01 binomial.

Four-RT result, re-run through the real analyzer → caller path:

| | Indura | Maxima | SSIV | TGIRT |
|---|---|---|---|---|
| Recall of any-channel signal site-obs, before → after | 38.7 → 100% | 85.4 → 100% | 71.3 → 100% | 58.8 → 100% |

The before figures are uneven partly because Maxima's deletions leaked into the old mismatch rate.
Sites with more than one row: 164 → 0. The 100% is partly circular, since the signal is defined by
the same thresholds; it checks wiring, not method.

## 5. Tests

- `tests/test_channel_priors.py` covers changes 1–3.
- `tests/test_channel_calling.py` covers change 4:
  - detection is the same whichever channel carries the signal
  - effect and significance gates both apply
  - the fitted background is calibrated on a simulated overdispersed null, where the binomial is not
  - priors never change the set of called sites
  - one row per site
  - `signature_type` overrides work
  - the analyzer's mismatch rate excludes deletions
- Full suite: 253 passed.

---

## 6. Site identity: why so many calls were `novel_candidate` (analysis before 4b)

The detection numbers above are sound. **The labels are not**, and there are three distinct causes. All
figures are four-RT, called sites, with MODOMICS relabelling enabled as the pipeline runs it.

**Label breakdown:**

| | Indura | Maxima | SSIV | TGIRT |
|---|---|---|---|---|
| labelled from profiles (`known`) | 5.0% | 17.2% | 19.0% | 16.5% |
| relabelled from MODOMICS (`known_modomics`) | 77.7% | 57.8% | 50.6% | 64.3% |
| `novel_candidate` | 17.4% | 25.0% | 30.4% | 19.1% |

The MODOMICS labels are acp3U 207, D 162, m1G 149, m7G 75, ms2i6A 42, I 34. Without MODOMICS
(an earlier harness run), about 84% of sites were novel. The pipeline passes MODOMICS, so these
numbers are the relevant ones.

The MODOMICS path is sound. `get_known_mods_linear` aligns the isodecoder's MODOMICS sequence to the
reference: 44 of 47 references have a map, and 99.7% of mapped modifications land on their parent
base. This is a sequence-specific, linear-coordinate map, with no Sprinzl arithmetic involved.

### Cause 1: hard-coded Sprinzl positions compared with linear positions

`ModificationProfile.typical_positions` holds canonical Sprinzl numbers, e.g. m7G `[46]` and m1A
`[58, 14, 9]`. The caller compares them directly with linear reference positions. Measured on the
MODOMICS map, the canonical position and the linear position differ like this:

| mod | linear − Sprinzl shift (count of tRNAs) |
|---|---|
| m7G46 | 0 (12), +1 (11), −1 (1) |
| m1G37 | +1 (6) |
| Ψ55 | −2 (2), −1 (4), 0 (13), +1 (12) |
| m5U54 | 0 (13), ±1–2 (18), **+9 to +18** (13; class II, long variable arm) |
| s4U8 | 0 (22), +1 (3), +2 (1) |

The shift varies by tRNA (D-loop insertions 17a/20a/20b, variable-loop length), so no single offset
fixes it. The effect on the 150 profile-labelled sites:

- **41 agree with MODOMICS**: Ψ 26/36, m7G 15/19.
- **All 68 `m5C` labels are wrong.** E. coli tRNA has no m5C. These are RT stops at linear 48, one
  base after acp3U47, which happens to match m5C's typical position 48.
- These wrong labels also **block** MODOMICS relabelling, which only applies to `novel_candidate`
  rows.

### Cause 2: RT stops are compared at the wrong base

A stop is recorded at the base 3′ of the modified one. **119 of 240 novel sites sit exactly one base
3′ of a MODOMICS modification**, almost all RT-stop-dominant. The caller compares the stop position
itself with the modification position, so it misses by one every time.

Offsets of called sites from the nearest known modification, by dominant channel:

| offset (mod − called) | deletion | mismatch | RT stop |
|---|---|---|---|
| 0 | 157 | 306 | 255 |
| −1 | 35 | 0 | **150** |
| −2 | 12 | 0 | 0 |
| none within 2 nt | 53 | 35 | 50 |

Substitutions land on the modified base. Deletions mostly do too, with a 3′ tail.

### Cause 3: genuinely unexplained

About 106 site-observations (10%) have no known modification within 2 nt. These are the real novel
candidates, or gaps in MODOMICS. The 3 references without a MODOMICS map also fall here.

### Not a cause: the Sprinzl `'20a'` crash

The standalone `modifications` command crashes in `annotate_signatures` (HEAD too). It calls that
function without `ref_seq`, which drops to the "merge on raw Sprinzl positions" branch and casts
`'20a'` to int. Stage 6 never takes that branch.

## 7. Implemented (change 4b)

1. **Signals are mapped to the modified base** (`_IMPLICATED_OFFSETS`), then identity is looked up
   there. Measured where only one of the two candidate bases is modified:
   - an RT stop at *p* implicates *p − 1* (58 cases) over *p* (11)
   - a deletion implicates *p* (17) over *p − 1* (6)
   - a substitution implicates *p* (60 to 1)

   The base used is reported as `modified_position`. `position` stays where the signal was
   observed.
2. **MODOMICS first.** A MODOMICS modification at an implicated base labels the site
   (`source='known_modomics'`, `identity_support='modomics'`). It overrides any profile label.
   Previously MODOMICS only relabelled novel rows. Every MODOMICS hit at the implicated bases leads
   `candidates`, followed by compatible profiles.
3. **Sprinzl `typical_positions` are no longer compared with linear positions.** A profile labels
   a site only on a specific substitution-pattern match (`identity_support='pattern'`). Otherwise
   the site is `novel_candidate`. `use_position_priors=False` disables the MODOMICS lookup.
4. **Offline MODOMICS for human and mouse.** `data/human_modomics_sequences.json` (42 sequences) and
   `data/mouse_modomics_sequences.json` (12) were exported from the MODOMICS cache. Stage 6
   defaults to `no_modomics: true`, so previously only E. coli had a map offline.

   Coverage with these files:

   | reference set | exact isodecoder | same-isotype donor | none |
   |---|---|---|---|
   | human | 294 | 151 | 12 |
   | RVHMS04_5215 library | 3,337 | 1,661 | 218 |
   | mouse | 50 | 40 | 151 |
5. **`'20a'` crash fixed.** The standalone CLI passes `ref_seq` and builds the same per-tRNA MODOMICS
   map as stage 6. `annotate_signatures` skips non-integer Sprinzl labels in its no-mapping
   branch instead of casting them.

Four-RT result (calls unchanged at 282/268/237/272):

| | Indura | Maxima | SSIV | TGIRT |
|---|---|---|---|---|
| `novel_candidate`, before → after | 17.4 → 12.8% | 25.0 → 15.3% | 30.4 → 16.9% | 19.1 → 13.6% |
| MODOMICS-labelled | 77.7 → 87.2% | 57.8 → 84.7% | 50.6 → 83.1% | 64.3 → 85.7% |

- **m5C labels fell from 68 to 0.**
- **26 former "Ψ" labels became m1G.** They are RT stops at 38 after m1G37; Ψ38 is RT-silent, and it
  stays in `candidates`.
- **Remaining novel sites: 154 site-observations** (deletion 68, RT stop 51, mismatch 35). These
  have no MODOMICS modification at the implicated base.

On the human example project, 48% of sites are MODOMICS-labelled.

### Open items from 4b

- **`unknown(O)` labels.** The MODOMICS symbol `O` is missing from the name table in `modomics.py`,
  giving 89–117 such labels on the human example.
- **Mouse coverage is thin** (151 of 241 references have no map). Human sequences could serve as
  cross-species donors.
- **`organism` defaults silently to E. coli.** Stage 6 uses `'Escherichia coli'` when the key is
  absent. The example project's config lacks it, so a human dataset was annotated with E. coli
  MODOMICS until `organism: human` was set.
- **Profile catalogue gaps** (acp3U, D, m5U, m2A, k2C) matter only for references with no MODOMICS
  map.

## 8. Novel sites after the MODOMICS lookup fixes

**Lookup fixes (`6a07e75`).** On the four-RT run, 129 of the 154 novel observations left after 4b
were known modifications the lookup had missed:

- **Isotype names.** `Ile2` and `fMet`/`iMet` didn't match MODOMICS keys. They now map to `ile`
  and `ini`.
- **Borrowed maps dropped the anticodon loop**, losing m1G37 in every E. coli tRNA-Pro. Stage 6
  now keeps these modifications, tagged `modomics_isotype`.
- **Wrong entries in the symbol table.** `?` is m5C (on C), not a modified G; `}` is k2C (on C), not
  a modified U. The table is now filled from MODOMICS's own symbol list (shipped as
  `data/modomics_symbols.json`).
- **Mitochondrial references** no longer borrow cytosolic maps.

Result: novel observations 13–17% → 3–7% by enzyme; calls unchanged.

**Novel-site table (`8e405c7`).** Stage 6 writes `results/modifications/novel_sites.{parquet,csv}`,
one row per unexplained site across samples. It records:

- channels fired and substitution spectrum
- median and maximum rate per channel
- the samples, conditions and enzymes showing the site
- labels the same site got in other samples
- the nearest MODOMICS modification and its offset

The modification report has a matching **Novel Sites** panel. On the four-RT run, 13 sites remain:

- **Within 3 nt of a known modification:** deletions 2–3 nt from ms2i6A37 / m1G38 on Maxima and
  SSIV, i.e. the signal smears beyond the offsets used for identity.
- **Real MODOMICS gaps:** e.g. `Pro-TGG`:34, probably cmo5U34, absent from the borrowed Pro-CGG
  map.
- **Single TGIRT RT stops** at about 120× coverage, probably noise.

**N-masked references: fixed with option (b).** `RTSignatureAnalyzer` computed the substitution
rate against an `N` reference base, so every read counted as a mismatch. On the human example
(`tRNA_database_masked`, 533 masked bases), 1,754 of 3,115 calls sat on masked bases at
mismatch ≈ 1.0. HEAD does the same.

**The signal at masked bases is mostly real.** The masked positions are 57% A and 37% G, and 242
of 512 coincide with a MODOMICS modification: m1A 89, m22G 38, I 25, m1I 20, m1G 17. These are
misincorporating modifications, masked because they disrupt mapping. Unmasked positions read
99.9% correctly. What was wrong was the measurement, not the detection.

**Fix:**

- **Restore the bases.** A new `unmasked_reference` config key (CLI `--unmasked-reference`) names
  the unmasked FASTA. `positional.unmask_reference` checks it has the same names and lengths and
  agrees at every unmasked base, then restores the masked bases before any rate is computed.
- **Fallback.** Without an unmasked reference, substitutions are not measured at N (deletion and
  RT stop still are), and stage 6 warns.

**Result on the human example:**

- All 533 bases restored; calls 3,115 → 3,113.
- Masked-base calls are now at their real rate (median 0.87).
- 778 of them are labelled: m1A 307, m22G 128, m1I 77, I 71, m1G 69.
- The 898 still novel are dominated by A→G and A→T, the inosine/m1I and m1A signatures, in tRNAs
  that MODOMICS's 42 human sequences don't cover.

## 9. Novel-site prediction (prototype, not in the pipeline)

`docs/prototypes/novel_site_classifier.py` trains a random forest on the exact-MODOMICS-labelled
sites and predicts labels for novel sites. Features are the site's RT fingerprint (median rate per
channel per enzyme × temperature), substitution spectrum and local sequence. Validation holds out
whole isotypes. Four-RT run, 88 labelled sites, 5 classes (D, ms2i6A, m7G, acp3U, m1G):

| features | held-out accuracy |
|---|---|
| majority class | 0.38 |
| reference base only | 0.73 |
| RT fingerprint + sequence | 0.84 |
| + anticodon-relative position | 0.88 |

Before the lookup fixes, it predicted m1G for `Pro-GGG`/`Pro-TGG`:38 (p 0.82–0.84). That is now
confirmed by MODOMICS. It also predicted acp3U for `Ile2`:47/48 and m7G for `Ile2`:46. It cannot
name unseen classes: `Ile2`:34 is k2C and got m1G at p 0.31. **Low probability must read as
"unknown".**

Why it stays a prototype:

1. **The fingerprint uses 12 RT × temperature conditions.** A normal run uses one enzyme;
   single-condition accuracy is unmeasured.
2. **Labels are only as good as MODOMICS.** There is no ground truth for the sites that matter (the
   unlabelled ones) and no true negatives.
3. **Small and narrow:** 88 sites, five classes, one organism.

What would answer these:

- **RT_comp (existing data) is enough to measure point 1.** Retrain on each enzyme × temperature
  subset alone and compare with the 12-condition model. No new experiment is needed.
- **A dedicated experiment would answer points 2 and 3.** The design would be E. coli wild type
  versus Keio single-gene knockouts of tRNA modification enzymes:
  - ΔtrmB (m7G46), ΔtapT (acp3U47), ΔmiaA/ΔmiaB (i6A/ms2i6A37)
  - ΔdusA/B/C (D), ΔtruA (Ψ38–40), ΔtrmA (m5U54)
  - an unmodified in-vitro-transcribed tRNA pool as the per-channel null

  Run them on two RTs with complementary channels (e.g. TGIRT and Maxima at 55 °C) in triplicate.
  Loss of signal in a knockout gives per-site ground truth and true negatives. The IVT pool also
  calibrates the background model directly. Essential enzymes (trmD for m1G37, tilS for k2C34)
  would need depletion strains.

## 10. Change 5: per-channel aggregates, grouped by site

`ReplicateAggregator` now groups calls by **site** (tRNA, position) within each condition, instead
of by (site, label). `aggregated_modifications` and `consensus_modifications` gain:

- **Per-channel means:** `mean_gap_rate` and `mean_rt_stop_pct` alongside `mean_mismatch_rate`.
- **Pattern and channel:** `dominant_pattern` and `dominant_channel` (mode over replicates), plus
  `channels_fired` (counts).
- **Labels:**
  - `labels` lists every label replicates gave, with counts.
  - `modification` is the most common resolved label; it is `novel_candidate` only if no replicate
    resolved one.
  - `source`, `identity_support`, `modified_position` and `rt_enzyme` describe the chosen label.
- **One p-value per replicate** enters the Fisher combination.

On the four-RT run, the 5 site×condition groups previously split by label are now single rows.
All five are acp3U47 vs m7G46: an RT stop at 47 points back to m7G46, a deletion stays on acp3U47.
3 of the 5 now reach consensus; the 2+1 split had denied it.

## 11. Change 6: regression test on real four-RT data

`tests/data/four_rt_fixture.parquet` (156 KB) holds per-position counts for 5 tRNAs × 36 four-RT
libraries, taken from stage 6's `mismatch_profile` + `rt_profile`:

- **Arg-ACG:** all channels.
- **Lys-TTT:** the spec's example, position 47.
- **Phe-GAA:** RT-stop dominant.
- **Pro-TGG:** deletion dominant.
- **Leu-GAG:** no deletions.

`tests/test_four_rt_regression.py` runs the stage 6 calling path on it in about 7 s. Following
decision D, every assertion is enzyme-agnostic:

| test | asserts |
|---|---|
| fixture exercises the problem | mismatch-only recall is uneven across enzymes (0.26–0.49) |
| recall even across enzymes | every enzyme ≥ 0.95, spread ≤ 0.05, and above mismatch-only |
| no regression | every substitution site is still called, with `mismatch` in `channels_fired` |
| calibrated on real null | fitted background: median ≤ 2% at p<0.01; single binomial rate: more than 2× that |
| enzyme label never changes detection | same sites with and without `rt_enzyme` |

**Checks on the test itself:**

- **It discriminates.** With deletion and RT stop disabled (a mismatch-only caller), recall is
  0.26–0.49 and the test fails.
- **It holds on the full dataset.** With `TRNASEQ_FOUR_RT_DIR` set, the recall test repeats on all
  49 tRNAs; it passed in 30 s. The spec's ordering ("Indura most, TGIRT least") is not asserted.
  It is a property of the baseline, not of the fix, and the Maxima/TGIRT gap was 3% (§3).

## CHANNEL_FIX_SPEC "Done when"

- [x] `rt_enzyme` / `rt_temp` flow from config to caller, `null`-safe (`824f1fa`)
- [x] `combined` is the default; per-modification override retained (`824f1fa`)
- [x] Per-enzyme priors derived from the four-RT data, not hand-set (`824f1fa`)
- [x] `combined` aggregation semantics documented, with a per-channel background model (§1–§5)
- [x] All three channel rates present in the aggregate outputs (`1597ee4`)
- [x] Regression test passes. The spec's enzyme ordering was replaced by enzyme-agnostic
  assertions (decision D).
- [x] Existing tests still pass (298 passed, 1 skipped)

Remaining open items are listed under §7 "Open items from 4b" and in §9.
