# Design: channel-aware calling and site identity (CHANNEL_FIX_SPEC change 4)

**Status (2026-10-03):** detection is decided and implemented (§1–§5). Site identity and position
handling (§6) is analysed with a fix proposed in §7; **it needs a decision before it is built.**
Changes 5–6 are not started.

Evidence comes from the 2024-06-19 four-RT run (4 enzymes × 3 temperatures × 3 replicates,
`ecoli.fa`, 49 references). Only `mismatch_profile.parquet` and `rt_profile.parquet` were used:
110,114 site-observations at coverage ≥ 100, with position 1 excluded.

Decisions taken:

| | Question | Decision |
|---|---|---|
| A | How channels combine | **OR** across substitution, deletion and RT stop |
| B | Role of per-enzyme channel priors | **Confidence and QC flag only**, never detection |
| C | Output shape | **One row per site** with a ranked candidate list |
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

## 6. Site identity: why so many calls are `novel_candidate`

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

## 7. Proposed fix (change 4b, not built)

1. **Map RT-stop evidence to the modified base.** For identity lookups, an RT stop at linear *p*
   implicates *p − 1*. Deletions use *p* (optionally *p + 1*). Detection is unchanged.
2. **MODOMICS first.** Look up the sequence-specific MODOMICS map at the implicated base before any
   profile ranking. A MODOMICS modification there:
   - becomes the label (`identity_support = 'modomics'`), overriding profile labels instead of only
     relabelling novel rows, and
   - is placed first in `candidates`.

   On the four-RT data, that combined with item 1 resolves about 119 of 240 novel sites and corrects
   the 68 false m5C labels.
3. **Real Sprinzl numbering for profile positions.** This is only needed where MODOMICS has no map,
   e.g. the 3 E. coli references, much of human/mouse, and the RVHMS04 library. Build a per-reference
   linear→Sprinzl table and compare `typical_positions` in Sprinzl space. Options:
   - **(a) Infernal/tRNAscan-SE, precomputed (recommended).** Run `cmalign` against the tRNA covariance
     model once per reference database, at `build-reference` time. Ship the table next to the FASTA;
     stage 6 reads it. This is standard structural numbering and handles D-loop and variable-arm
     insertions. It is installed locally (`/usr/bin/cmalign`, `/usr/local/bin/tRNAscan-SE`).
     **Whether it is on O2 is unknown.** Precomputing means O2 only needs the table.
   - **(b) Pure-Python anchored numbering.** Use three anchors:
     - 5′ end for 1–9
     - anticodon (34) for 27–43
     - 3′ CCA (76) for 49–76

     The D-loop and the 44–48 variable region stay ambiguous and would be flagged. This has no
     dependency but is approximate exactly where m7G46/acp3U47 sit.
4. **Fix the `'20a'` crash.** The CLI passes `ref_seq` to `annotate_signatures`, as stage 6 does. The
   raw-Sprinzl merge skips non-integer positions instead of crashing.
5. **Profile catalogue (optional, later).** acp3U, D, m5U, m2A and k2C have no profile. That matters
   only where MODOMICS has no map.

Expected outcome on four-RT: about 10% `novel_candidate` (cause 3 only), and no profile label
contradicting MODOMICS.

### Decisions needed

- **E. Scope:** build change 4b (items 1, 2, 4) now, before change 5?
- **F. Sprinzl numbering:** (a) precomputed with Infernal, or (b) pure-Python anchors? Is Infernal
  available on O2, or is precompute-locally-and-ship acceptable?
- **G.** CHANNEL_FIX_SPEC said not to touch Sprinzl handling in this change. 4b does touch it, in the
  identity layer only; detection stays unchanged. Should it stay a separate commit (recommended), so
  coordinate changes can be reverted on their own?

---

Not yet done: change 5 (per-channel rates and `dominant_pattern` in the aggregates; also make
`ReplicateAggregator` group by site rather than by label, since 5 of 447 site×condition groups split
when replicates disagree) and change 6 (enzyme-agnostic regression tests on a four-RT fixture).
