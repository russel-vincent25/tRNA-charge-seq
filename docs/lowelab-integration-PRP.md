# PRP: Lowe Lab tRNA tool integration

**Status:** proposal / backlog. Nothing here is implemented.
**Scope of this document:** planning only. No code changes accompany it.

## How to read this

Every claim below cites `file:line` against the repo as of this commit. Line
numbers were spot-checked when the document was written; if code has moved,
re-verify before acting.

Items marked **[verify]** are readings of the source that have *not* been
confirmed against real data or a live run. Treat them as leads, not bug
reports. Do not "fix" a **[verify]** item before reproducing it.

---

## 1. Positioning: why integrate at all

The load-bearing observation:

> **No Lowe Lab tool measures aminoacylation.**

tRNAscan-SE, GtRNAdb, tRNAviz, tDRnamer and tRAX are about *finding,
naming and comparing* tRNAs. Charge measurement — the 3' CA/CC
discrimination this pipeline is built around (`src/stats_collection.py:195-198`,
`src/plotting.py:204-207`) — is absent from all of them. tRAX, the closest
counterpart pipeline, does alignment, coverage, fragment classes,
misincorporation and differential expression, but not charging.

So the two toolsets are complementary. Charge is this repo's contribution;
theirs is **annotation and standardisation**. The integration thesis is: adopt
their vocabulary and coordinates so this repo's charge data becomes
comparable and citable, without giving up anything it does uniquely.

| Lowe Lab tool | What it would give us | State in this repo today |
|---|---|---|
| **tRNAscan-SE** | tRNA gene finding → the reference DB | Never invoked. Only its *precomputed output* is consumed, and only for intron removal (`tRNA_database/tRNAdb_tool.py:34-61`) |
| **GtRNAdb** | Mature sequences + canonical IDs | Already the source; our ID grammar is GtRNAdb-derived |
| **tDRnamer** | Standard tDR naming ontology (`tDR` + Sprinzl range + GtRNAdb ID) | Ad-hoc `is_5p_tRF` / `is_degraded_tRNA` / `is_pre_tRNA` (`src/stats_collection.py:170-178`) |
| **tRNAviz** | Sprinzl canonical position numbering | None anywhere in the codebase |
| **tRAX** | Read-specificity coverage classes | Isotype / codon / transcript levels exist; **no isodecoder level** |

### Provenance note

`trna.ucsc.edu` and `lowelab.ucsc.edu` were **blocked by the network egress
proxy** in the environment where this document was written. Tool descriptions
were sourced from the `UCSC-LoweLab` GitHub organisation and from web search,
**not** from the canonical site. Anyone extending this document should
re-check the primary source, in particular the exact tDR naming grammar and
the covariance model filenames.

---

## 2. Backlog

Ordered so that each tier can be done without the ones below it.

### Tier 1 — Defects to resolve before building on top

These were found while surveying the code for integration points. They are
listed first because Tiers 2-4 touch the same code paths.

- [ ] **`align_3p_nt` vs `align_3p_nts` column mismatch. [verify]**
  `src/transcript_mutations.py:48` declares a column `align_3p_nt`, and
  `:276` filters on it comparing against *single* characters:
  ```python
  ((sample_stats['align_3p_nt'] == 'A') | (sample_stats['align_3p_nt'] == 'C'))
  ```
  But `src/stats_collection.py:39` writes `align_3p_nts` (plural), holding
  *two*-character values — `'CA'`, `'CC'`, `'GA'`, `'CG'`
  (`stats_collection.py:330`, `:195-198`). Every other consumer uses the
  plural name (`src/plotting.py:48`, `:204-207`, `:1036`).
  Reading suggests `find_muts()` would raise `KeyError` on stats CSVs
  produced by the current code, and that even with the name corrected the
  single-char comparison would never match a 2-char value.
  Looks like a regression from the enhanced-branch 3'-end change.
  **Action:** run `find_muts()` against a real `*_stats.csv.bz2` first and
  capture the actual failure before changing anything.

- [ ] **Mixed units in the fragment aggregation. [verify]**
  `src/stats_collection.py:212-216` aggregates the fragment flags by summing
  booleans:
  ```python
  Total_5p_tRFs=("is_5p_tRF", "sum"),
  ```
  Summing a boolean column counts **rows** — i.e. distinct read species,
  since each row carries its multiplicity in `count`. The 3'-end totals
  alongside them use `("count", "sum")` (`:200-203`), i.e. **reads**. Both
  are then merged into the same `agg_df` (`:206-209`, `:218`), so columns in
  one table are on two different scales.
  Second discrepancy: the fragment flags are computed on `stat_df` (all
  reads) while the 3'-end totals use `filtered_stat_df` (full-length only,
  `:185-186`) — different denominators.
  **Action:** decide whether these columns should be read counts (likely) and
  whether they should share a denominator with the charge columns.

- [ ] **Fragment categories are not mutually exclusive. [verify]**
  `is_pre_tRNA` (`src/stats_collection.py:178`) is an OR over non-templated
  ends, so it can be true simultaneously with `is_5p_tRF` or
  `is_degraded_tRNA`. Meanwhile `src/plotting.py:232-236` puts
  `5p_tRF + degraded` into the charge denominator:
  ```python
  100 * row['CA_count'] / (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded'])
  ```
  If those categories can overlap, the denominator double-counts.
  **Action:** confirm the intended exclusivity, then either make the classes
  disjoint or document the overlap explicitly. This matters for any published
  charge number, so it should be settled before Tier 3 renames these fields.

- [ ] **Single position index applied to all transcripts.**
  In `plot_transcript_mut_pos`, `idx_pos` is computed once at
  `src/transcript_mutations.py:1170`:
  ```python
  idx_pos = (tr_pos - idx_start)
  ```
  but the per-annotation loop only begins at `:1221`, and `mut_mat[:, idx_pos]`
  at `:1266` then uses that one scalar for every transcript. Combined with the
  Tier 2 problem below, the plotted position is not the same biological
  position across rows. **Action:** move the resolution inside the loop.

### Tier 2 — Sprinzl canonical numbering (the foundation)

**Problem.** There is no canonical tRNA numbering anywhere in the codebase.
Positions are raw indices into each transcript's own sequence.

`plot_transcript_mut_pos(tr_pos=34)` documents `tr_pos` as *"Transcript
position, counted from left to right"* (`src/transcript_mutations.py:1154`).
That is a 5'-end index, **not** Sprinzl 34. It coincides with the wobble base
only for transcripts with canonical D-loop and variable-loop lengths, and it
is shifted by one for every His transcript, which carries a prepended `G-1`
(`tRNA_database/tRNAdb_tool.py:82-83` for mito, `:110-111` for cyto). The notebooks call it
with `tr_pos=34` (`projects/tRNAseq_mtests/process_data.ipynb`) and
`tr_pos=37` (`projects/example/process_data.ipynb`), which suggests the
intent is structural even though the mechanism is positional.

The only cross-transcript devices that exist today are:

1. **Zero-padding right-aligned to CCA** — `src/transcript_mutations.py:602-605`,
   `:702-704`: *"Insert from right side so CCA is always indexed at the
   highest index"*. Correct for the acceptor stem and T-arm, wrong as soon as
   D-loop or variable-loop lengths differ. The CSV header written at
   `:457` is `'P'+str(i+1)` over `self.longest_tRNA`, so `P51` means "column
   51 of the longest-tRNA frame", not "nucleotide 51 of this tRNA".
2. **A proportional percentile stretch** — `src/plotting.py:1051-1062`:
   ```python
   len_map_len[len_i] = np.percentile(np.arange(max_len), np.linspace(0, 100, len_i), method='nearest')
   ```
   applied at `:1078`. This is a length-only rescale keyed on
   `tRNA_annotation_len`. It is exactly the slot a real Sprinzl table
   should occupy.

**Proposal.** Compute Sprinzl positions **once, at database-build time**, and
ship them as a sidecar table beside each reference FASTA.

Both tDRnamer and tRNAviz derive Sprinzl numbering the same way: Infernal
`cmalign` against a covariance model specialised for alignment and numbering.
`UCSC-LoweLab/tRNAviz-data` contains a reference implementation
(`tRNA_position.py`, driven by `parse-tRNAs.py -n <model>`) and emits a table
with one column per position.

Doing this offline at build time is the key design choice: it keeps Infernal
out of the analysis-time dependency set, and avoids depending on tDRnamer
itself, which requires **Python 2.7** plus bowtie2, BLAST+, EMBOSS, samtools
and Infernal — a hard fit for this Python 3 codebase.

- [ ] Generate a per-transcript `index → Sprinzl label` sidecar TSV per
      reference DB, alongside the existing `.phr/.pin/.psq` BLAST files.
- [ ] Inject it at `src/transcript_mutations.py:105-118` — the `SeqIO.parse`
      loop that builds `seq`, `seq_len`, `PSCM`, `RTstops`, `mut_freq`,
      `gap_freq`. Everything downstream reads
      `tr_muts_combi[species][anno][...]`, so one addition there reaches all
      consumers.
- [ ] Mirror it in `src/misc.py:31-44` (`read_tRNAdb_info`) so the stats and
      plotting sides can reach it too.
- [ ] Convert consumers, cheapest and highest-value first:
      `plot_transcript_mut_pos` (`:1170`, `:1266`, title `:1302`), then
      `write_transcript_mut` header (`:457`), `plot_transcript_cov`,
      `plot_transcript_mut`, `plot_transcript_mut_compare`,
      `plot_transcript_mut_cluster`, and `src/plotting.py:1051-1078`.

**Known blockers — resolve these as part of the work, not after:**

- `src/plotting.py:879-881`: `aa_cov_cols` contains
  `['sample_name_unique', 'tRNA_annotation_len', 'align_5p_idx',
  'align_3p_idx', 'AA_letter', <RPM col>]` — **no `tRNA_annotation`**. The
  frame is grouped by these at `:1049`, so transcript identity is aggregated
  away before any position mapping happens. A per-transcript Sprinzl lookup
  is impossible until `tRNA_annotation` is added to that list.
- The reference FASTA has **introns already removed** and His carries a
  prepended `G-1` (intron cut at `tRNA_database/tRNAdb_tool.py:104-109`; His at `:82-83` and `:110-111`). Any
  Sprinzl coordinate obtained from GtRNAdb or tRNAscan-SE must be
  intron-adjusted and His-adjusted to index the FASTA the pipeline actually
  aligns against.
- `src/transcript_mutations.py:1330` (`freq_mut[~min_count_mask] = 0`)
  mutates arrays **in place inside `tr_muts_combi`**. Any remapping must not
  alias these buffers.
- Hardcoded 3'-anchored assumptions to re-express in canonical terms:
  `_sort_anno` uses `PSCM.sum(1).values[-2]` (`:1346`, the `C` of `CCA` —
  Sprinzl 75) and the end-fix at `:368-374` overwrites row `seq_len-1`
  (Sprinzl 76).
- `self.longest_tRNA` (`:103-110`) currently sets plot axis width; it would
  become a fixed canonical axis (1-76 plus the `17a/20a/20b/e11-e17`
  insertions) rather than an observed maximum.
- `src/simulate.py:53-58` hardcodes a plainly Sprinzl-derived `mod_sites`
  list applied as **raw** indices. It is self-consistent only because
  `_gen_ref` filters references to exactly 76 nt. Same category error as
  `tr_pos=34`, quarantined by a length filter. Fold it into the canonical
  scheme.

### Tier 3 — tDRnamer-compatible fragment naming

**Depends on Tier 2** (the Sprinzl range is half the name).

Goal: replace `is_5p_tRF` / `is_degraded_tRNA` / `is_pre_tRNA` with the
published ontology — `tDR` prefix, Sprinzl range, GtRNAdb ID — so fragment
output is interoperable with the wider tDR literature.

**The hard constraint.** This pipeline is bound to a positional,
`-`-delimited transcript name whose prefix must contain **no hyphen**:

```
<PREFIX>-<AminoAcid>-<ANTICODON>[-<isodecoder>-<gene>]
   [0]        [1]         [2]
```

`src/misc.py:41-43` is the single derivation point:
```python
tRNA_data[record.id]['codon']      = anticodon2codon(record.id.split('-')[2])
tRNA_data[record.id]['anticodon']  = record.id.split('-')[2]
tRNA_data[record.id]['amino_acid'] = record.id.split('-')[1]
```
This is why `tRNA_database/tRNAdb_tool_ecoli.py:14` rewrites
`Escherichia_coli_str_K-12_substr` → `..._K12_substr`: one stray hyphen would
shift every field.

Additional fragilities that a hyphen-rich tDR name would hit:

- Closed isotype vocabulary: the `AAA2A` dict (`src/plotting.py:1126`, applied
  at `:113` and `:1031`) raises `KeyError` on an unknown amino-acid token.
- `src/plotting.py:111` unconditionally strips a trailing `1`/`2` from the
  amino-acid token — intended for `Leu1`/`Leu2`, destructive for anything
  else numbered.
- `anticodon2codon` (`src/misc.py:26-28`) is ATGC-only; `NNN` raises.
- `@` is reserved as the multi-mapping join separator
  (`src/alignment.py:408`), so it must never appear in an ID.
- Inconsistent substring sniffs for compartment and spike-ins:
  `'mito' in n` (`src/alignment.py:413`), `'mito_tRNA' in anno`
  (`src/plotting.py:151`), `'Synthetic' in anno` (`:152`),
  `'Escherichia_coli' in anno` (`:1028`).

- [ ] **Do not loosen the ID parse.** Introduce a sidecar annotation table
      that replaces `read_tRNAdb_info` (`src/misc.py:31-44`) as the single
      place where isotype / anticodon / codon are derived. Parsing richer
      names positionally will silently corrupt fields rather than fail.

### Tier 4 — Independent of Sprinzl

- [ ] **Invoke tRNAscan-SE directly to build the reference DB.**
      Today the DB build `wget`s precomputed files from a personal
      `mim-tRNAseq` fork (`tRNA_database/README.md:14-57`), which is fragile and
      limits the pipeline to the species someone has already staged.
      Worth noting while here: **only the human build passes the intron
      arguments** (`tRNA_database/README.md:25`). The mouse (`:36`) and yeast
      (`:47`) invocations omit both, and E. coli uses a separate script
      (`:57`) — so yeast introns are never removed by this pipeline.
- [ ] **tRAX-style read-specificity classes** — transcript / isodecoder /
      isotype / non-specific — so ambiguous multi-mapping reads are reported
      honestly instead of collapsed. Today `src/plotting.py:248-313` provides
      isotype (`aa`), codon and transcript levels, but there is **no
      isodecoder level**: the `-1-1` isodecoder/gene fields are never parsed
      anywhere in `src/`. Related: for a multi-mapped read only the first
      alphabetically-sorted hit sets the AA/codon
      (`src/stats_collection.py:318`, `:356-358`), policed after the fact by
      the `single_codon` / `single_aa` flags (`src/plotting.py:130-150`).
- [ ] **Make the spike-in control sequences configurable.**
      `tRNA_database/README.md` already flags this: *"Note that the control
      sequences are hard-coded in the `tRNAdb_tool.py` script..."*
      (`tRNA_database/README.md:11`). They live
      at `tRNA_database/tRNAdb_tool.py:15-19`. A `--controls <fasta>` /
      `--no-controls` pair defaulting to the current built-in E. coli
      Lys/Thr pair keeps existing output byte-identical, and lets
      `tRNA_database/Synthetic_control_sequences.fa` be used without editing
      the script. Small, self-contained, no dependencies.

---

## 3. Non-goals

- **Do not replace SWIPE with bowtie2** to match tRAX. The `--symtype 1`
  protein-mode trick plus a custom substitution matrix
  (`utils/nuc_score-matrix_2.txt`) is deliberate and documented in
  `tRNA_database/README.md`: SWIPE only accepts a custom matrix in protein
  mode, and that is what makes `N`-masking of modified positions possible.
  The whole masked-database strategy
  (`src/transcript_mutations.py:1406-1624`) depends on it. Adopting Lowe Lab
  *vocabulary* does not require adopting their *aligner*.
- **Do not vendor tDRnamer.** Python 2.7. Use the covariance-model approach
  offline instead (Tier 2).
- **Do not reimplement charge measurement** to match any Lowe Lab tool.
  None of them do it.

## 4. Open questions

1. Which reference genomes must the Sprinzl sidecar cover? The shipped DBs
   include human, mouse, yeast, E. coli plus several project-specific
   libraries (`RVHMS*`, `RusselLib`, `YLDC`) whose IDs use a
   `tdbD…_tRNA-…` shape. Do the custom libraries need canonical numbering,
   or only the model organisms?
2. Should the masked databases (`tRNA_database_masked/`) carry their own
   sidecar, or inherit the unmasked one? IDs and lengths are preserved by
   `write_masked_tRNA_database` (`src/transcript_mutations.py:1557-1600`),
   so inheriting looks safe — worth confirming.
3. Is the tDR ontology wanted in *output only*, or should it also become the
   internal fragment representation? Output-only is far cheaper and avoids
   the ID-grammar problem entirely.
