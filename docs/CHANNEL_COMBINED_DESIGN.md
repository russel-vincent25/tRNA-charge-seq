# Design: channel-aware calling and site identity (CHANNEL_FIX_SPEC change 4)

**Status (2026-10-03):** detection (change 4, §1–§5), site identity (change 4b, §6–§7) and the
novel-site lookup fixes and table (§8) are implemented. Prior art and two corrections from the
2026 RT literature are in §12. The novel-site predictor (§9) is a
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
RT stop 0.20 at the time of this analysis; now 0.10, see §12) **and** its BH q < 0.01. The
false-positive proxy is replicate
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
  SSIV. ⚠ **Superseded by §12:** this is homopolymer indel-placement ambiguity, not a signal smear —
  the modification lies inside the same homopolymer run as the called deletion.
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

## 12. Prior art: the channel model is published, and our offset rule is confirmed

Added 2026-10-03 after reading four papers (PDFs supplied by RV; metadata and DOIs via PubMed).
**The channel-aware model is not novel — it is the field's standard, and this pipeline was behind
it.** What follows is what each paper settles, and the two places our implementation disagreed.

| paper | channels scored | temperature | enzymes | organism | replicates |
|---|---|---|---|---|---|
| **Nakano et al. 2025**, *Nat Commun* 16:1047, [doi](https://doi.org/10.1038/s41467-025-56348-1) — the primary paper | misincorporation + RT stop | 25/37/42/55 °C x 1/2/16 h, **Induro only** | Induro vs TGIRT, **42 °C only** | human, mouse | n = 2 technical (temp series) |
| **Nakano, Gamper & Hou 2026**, *Methods Enzymol* 725:51, [doi](https://doi.org/10.1016/bs.mie.2025.10.007) — methods chapter summarising the above | same | same | 6 RTs for read-through, **each at its own published condition** | human | same |
| **Pedor et al. 2026**, *RNA Biol* 23:1, [doi](https://doi.org/10.1080/15476286.2026.2720028) | **mismatch only** | fixed 42 °C / 16 h | MRT, MRT-CBD, uMRT, Induro | yeast | n = 6 technical |
| **Werner et al. 2020**, *NAR* 48:3734, [doi](https://doi.org/10.1093/nar/gkaa113) | misincorporation, arrest, **nucleotide skipping** | — | 13 RTs | model RNAs | — |

### Confirmed: the +/-1 window and three per-modification behaviours

Nakano et al. state it with an exhaustive search:

> "At each of these modifications, defined as position 0, we only found RT stops or RT
> misincorporation from -1 to +1, despite an exhaustive search from positions -2 to +2, indicating
> no RT jump-over."

They also report that **acp3U20a and ms2i6A37/ms2t6A37 "responded primarily by RT stops... at the
+1 position after the site of the modification"**, while **I34 responded "exclusively by
misincorporation without RT stop"**. All three match what change 4b derived independently from the
four-RT data, and `_IMPLICATED_OFFSETS` encodes exactly that window. **Cite Nakano for the offset
rule; do not present it as new.**

Verified per channel on our 1,059 calls (offset = known modification position - observed signal
position):

| offset | deletion | mismatch | rt_stop |
|---|---|---|---|
| -3 | 0 | 0 | 2 |
| -2 | 12 | 0 | 0 |
| **-1** | 35 | 0 | **173** |
| **0** | **192** | **341** | **266** |
| +1 | 5 | 1 | 5 |
| +2 | 4 | 0 | 0 |
| +3 | 12 | 0 | 0 |

- **Mismatch: zero calls outside +/-1.** Exact agreement.
- **RT stop:** the 173-call pile at -1 *is* their "stop at the +1 position after the modification".
  Six sites fall outside: four are single observations at 108-155x on TGIRT only (noise), one is
  marginal (`Ser-TGA`:41), and one is real — see below.
- **Deletion:** three sites appeared to violate the window. **All three resolve, and the earlier
  claim in this document that deletion signal "smears 2-3 nt" was wrong.**

### Retracted: deletions do not smear. It is homopolymer indel placement.

| site | apparent offset | cause |
|---|---|---|
| `Pro-TGG`:36 | +2 to m1G38 | **m1G38 lies inside the same GGG run (36-38)** — a deletion within a homopolymer run has no unique alignment |
| `Trp-CCA`:39 | -2 to ms2i6A37 | **ms2i6A37 lies inside the same AAAA run (36-39)** |
| `Pro-TGG`:35 | +3 to m1G38 | linear 35 **is** the wobble base (Sprinzl 34): unannotated cmo5U34 |

Once homopolymer ambiguity and that one missing annotation are accounted for, **Nakano's +/-1
window holds in our data across all three channels.**

✅ **Implemented 2026-10-03.** `_homopolymer_run()` and `_implicated_positions()` replace the bare
offset table for the deletion channel: a deletion at *p* implicates *p*, *p-1*, and then every base
of the homopolymer run containing *p*, nearest first. Sites where the run is longer than one base
carry **`position_ambiguous=True`**, so an ambiguous placement is visible rather than silent. The
substitution and RT-stop channels are untouched, and with no reference sequence the behaviour falls
back to the offsets.

On the four-RT data, **detection is identical (1,417 call-observations either way)** and 19 calls
across exactly the three predicted sites are relabelled:

| site | was | now | modified_position |
|---|---|---|---|
| `Trp-CCA`:39 (12 obs) | novel_candidate | **ms2i6A** | 37 |
| `Pro-TGG`:36 (4 obs) | novel_candidate | **m1G** | 38 |
| `Arg-ACG`:46 (3 obs) | novel_candidate | **m7G** | 47 |

Unique sites 96 -> 95 (two merge onto one modified base), MODOMICS-explained 95.1% -> 96.5%, and
novel_candidate 4.9% -> 3.5%. 227 of 1,417 calls are flagged `position_ambiguous`, 159 of them
deletion-dominant.

### The one genuinely unexplained reproducible site

`Pro-TGG`:34 — RT stop 58% median, **7 observations across Indura and TGIRT at 2,275x**, no
MODOMICS modification within 3 nt. It sits in the anticodon loop of a tRNA whose map is borrowed
from Pro-CGG and is already known to be missing cmo5U34. Most likely unannotated anticodon-loop
chemistry; a *cmoB* knockout would settle it.

### Threshold: ours is stricter than the published standard, and it costs sites

Nakano et al. call a modification at **">10% of RT misincorporation or stop"**, with unmodified
positions 70/74/75 as controls. Our defaults are mismatch 0.10 (matching), deletion 0.10, and
**RT stop 0.20 — twice theirs**. Re-running the four-RT data at both:

| rt_stop threshold | calls | unique sites | MODOMICS-explained | rt_stop-dominant (explained) | sites in >=2/3 reps |
|---|---|---|---|---|---|
| 0.20 (old) | 1,059 | 82 | 95.1% | 457 (96.1%) | 81.7% |
| **0.10 (adopted)** | **1,417** | **96** | **95.1%** | 836 (95.8%) | 81.8% |

⚠ **Correction.** This row first read 1,411 calls / 830 stop-dominant. The sensitivity script
overrode the *calling* threshold but left `estimate_channel_backgrounds` fitting its null bulk at
`rate < 0.5 x` the **old** threshold, so the background and the gate disagreed. With both
consistent at 0.10 the figure is 1,417. The direction and size of the effect are unchanged.

**+33% calls and +14 sites at identical precision** on both proxies. The 14 gained sites are
exactly the classes we were under-detecting, all RT-stop-driven: **D16 in `Asn-GTT` (5 obs, 3
enzymes) and `Val-GAC-2` (9 obs, 2 enzymes), m7G46 in `Thr-GGT` (4 obs, 2 enzymes), s2C33, and a
Psi40** — dihydrouridine was our worst class (23/83) and Psi was 0/76. The other 9 gained sites are
single observations and look like noise.

#### Adopted 2026-10-03: `DERIVATION_THRESHOLDS['rt_stop']` 0.20 -> 0.10

CHANNEL_FIX_SPEC's "do not lower thresholds to compensate" was about substituting stringency for a
channel fix. This is the opposite case: after the channel fix, aligning with the published standard,
with no measured precision loss. Consequences, all verified:

- **Priors regenerated** from the same four-RT profiles. Signal observations rise (Indura 282 ->
  398, Maxima 268 -> 331, SSIV 237 -> 309, TGIRT 272 -> 379).
- **RT stops now lead for three of four enzymes** pooled (Indura 0.74, SSIV 0.49, TGIRT 0.54;
  Maxima deletion 0.35). A lower stop threshold admits many low-level stops, so "which channel tops
  out" is no longer the informative summary. **The enzyme contrast survives as relative
  differences**, which is what the fix was about: Indura has the highest stop weight (0.74) and the
  lowest deletion weight (0.07), Maxima the highest deletion weight (0.35), and Maxima at 55 °C
  remains the deletion-dominant cell of the matrix (0.45 deletion vs 0.16 stop).
- **`unexpected_channel` now fires on different cases.** Maxima at 55 °C no longer trips it (stop
  weight 0.036 -> 0.162, above the 0.10 cut-off); a deletion-only site on Indura does
  (weight 0.067), which is the biologically sensible flag — Indura barely deletes.
- **Detection ceiling barely moves: 65 -> 70 of 343** known modified bases (19% -> 20%), gaining
  D +2, Psi +1, m7G +1, s2C +1. The modifications still missed (Psi 75/76, m5U 44/44, s4U 28/28,
  t6A, m2A, cmo5U, Q) leave **no** RT signal rather than a sub-threshold one, so no threshold
  recovers them.
- **Three tests were rewritten**, not weakened: the priors test now asserts the relative enzyme
  contrast rather than which channel tops out; the `unexpected_channel` test uses Indura-deletion
  instead of Maxima-stop; and a pre-existing fold-change test now fires a single channel so the
  dominant channel is unambiguous. The four-RT regression test was re-verified to still
  discriminate — honest run 1.00 recall for all four enzymes, mismatch-only mutation 0.23-0.42
  and failing.

### What RT_comp still adds

Stated conservatively, since three of the four papers post-date the experiment:

1. **A factorial design.** No paper crosses enzyme with temperature: Nakano varies temperature for
   Induro alone and compares Induro vs TGIRT at 42 °C only; the Hou chapter runs each of 6 RTs at
   its own published condition, so enzyme and temperature are confounded there; Pedor fixes 42 °C.
2. **Biological triplicates**, against n = 2 technical for the temperature series.
3. **The deletion channel.** Absent from all three tRNA-seq papers (in Nakano the only "indel"
   mentions are `--no-indels` in cutadapt adapter trimming). **But deletions as a signature are
   established in the wider field** — Werner et al. list "nucleotide skipping" as one of three
   readouts, and BID-seq deliberately converts Psi into deletions. The defensible claim is narrow:
   *deletion as a major, enzyme- and temperature-dependent channel in tRNA-seq* (Maxima at 55 °C
   carries ~51% of its signal weight there, and a misincorporation+stop framework misses it).
4. **Bacterial tRNA** — the others are human, mouse and yeast.

Their detection ceiling matches ours: they cannot detect Psi, and reach m7G only with a lowered
4-12% cut-off (we found Psi 0/76, m7G 12/24). They note Marathon reaches both at 0.1-1%
misincorporation, so that ceiling is enzyme-specific, not fundamental.

### Precedent for the classifier prototype (§9)

Nakano et al. **cite Werner et al.** (their ref 22) but build no model. Their approach to ambiguous
sites is manual cross-reference of two RT datasets:

> "having two datasets of tRNA modifications, each with a different RT in a different workflow,
> would help resolve ambiguity, strengthen the prediction, and provide the basis for cross-reference
> of each dataset."

That is the same idea as the multi-condition fingerprint in §9, done by hand with two enzymes
instead of learned over twelve conditions. Werner et al. is the direct precedent for the learned
version. Both belong in the prototype's docstring.

### Mechanism worth keeping

> "the average misincorporation rate through all tRNA sequences remained constant at ~3%...
> indicating that the increase in readthrough was driven by decreases of RT stops."

Read-through gains come from losing RT stops while misincorporation holds constant. Channel
partitioning is therefore a function of condition, not only of enzyme — which is why the priors are
keyed on enzyme *and* temperature, and why they must never gate detection.

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
