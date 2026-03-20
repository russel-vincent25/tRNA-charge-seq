# tRNA Modification Analysis Module

Detects and quantifies tRNA modifications from RT (reverse transcription) signature analysis. Integrates with the MODOMICS database for annotation of known modifications and supports discovery of novel modification sites.

## Overview

The module identifies positions in tRNA transcripts that show elevated mismatch rates, gaps, or RT stops — all indicators of post-transcriptional modifications. It runs as Stage 6 of the trnaseq pipeline or as standalone CLI/Python commands.

### Pipeline Integration (Stage 6)

When `run_modification_analysis: true` is set in the config, Stage 6 runs 7 phases:

1. **PSCM extraction** — Builds Position-Specific Count Matrices from SWalign JSONs (parallelized with mpire)
2. **Background estimation** — Estimates per-sample error rates from low-modification positions
3. **Modification calling** — Tests each position against background using binomial test + Benjamini-Hochberg FDR
4. **Replicate aggregation** — Merges calls across replicates (majority vote + Fisher meta-analysis)
5. **Summary** — Generates per-tRNA modification summary tables
6. **SLAC crosstalk** — Analyzes per-read modification coordination using Fisher's exact test (parallelized)
7. **Report** — Interactive HTML dashboard with modification profiles

## Module Components

```
trnaseq/modifications/
├── positional.py          # PositionalExtractor — PSCM from SWalign JSONs
├── rt_signatures.py       # RTSignatureAnalyzer — mismatch/gap/RT-stop analysis
├── modification_caller.py # ModificationCaller — statistical calling + profiles
├── modomics.py            # MODOMICSAnnotator — MODOMICS API + fallback CSVs
├── crosstalk.py           # CrosstalkAnalyzer — SLAC per-read coordination
├── reference_builder.py   # MODOMICS isodecoder reference FASTA generator
└── data/                  # Shipped fallback data
    ├── ecoli_known_modifications.csv
    ├── human_known_modifications.csv
    ├── mouse_known_modifications.csv
    └── ecoli_modomics_sequences.json
```

### PositionalExtractor (`positional.py`)

Streams SWalign JSON files and builds a Position-Specific Count Matrix (PSCM) for each tRNA:

- Matrix dimensions: `ref_length × 8` (A, C, G, T, N, gap, coverage, rt_stop)
- Supports parallel processing across samples via mpire
- Handles bz2-compressed JSON files
- Filters by alignment quality (Fmax_score threshold)

### MODOMICSAnnotator (`modomics.py`)

Integrates with the MODOMICS database to annotate known modification positions:

- **API fetch**: Queries MODOMICS REST API for organism-specific tRNA modifications
- **JSON cache**: Caches API results locally (v2 format with raw sequences for alignment)
- **Shipped fallback**: Uses bundled CSVs + JSON when API is unavailable
- **Alignment-based mapping**: Uses `Bio.Align.PairwiseAligner` to map MODOMICS positions to user's reference FASTA (handles sequence differences between MODOMICS and custom references)
- **Isotype-level transfer**: When MODOMICS has no entry for a specific tRNA (isotype, anticodon), transfers modifications from a same-isotype donor, excluding anticodon loop positions (Sprinzl 32–38)

### ModificationCaller (`modification_caller.py`)

Statistical calling of modifications with 23 supported profiles:

- Binomial test against per-sample background error rate
- Benjamini-Hochberg FDR correction (returns boolean significance array)
- Position-specific priors from MODOMICS annotations
- Replicate aggregation with `ReplicateAggregator` (majority vote + Fisher meta-analysis)

### CrosstalkAnalyzer (`crosstalk.py`)

Single-read Analysis of Crosstalks (SLAC) following Behrens et al. 2023 (NAR):

- For each pair of known modification positions on a tRNA, builds a 2×2 contingency table from individual reads
- Fisher's exact test yields odds ratio and p-value
- OR > 1: modifications co-occur (positive coordination)
- OR < 1: modifications are anti-correlated (negative coordination)
- Parallelized across samples with mpire

## Supported Modifications

| Modification | Full Name | RT Signature |
|--------------|-----------|--------------|
| **m1A** | 1-methyladenosine | A→any mismatches + strong RT stops |
| **m1G** | 1-methylguanosine | G→any mismatches + RT stops |
| **m3C** | 3-methylcytosine | C→T mismatches + RT stops |
| **m5C** | 5-methylcytosine | C→T mismatches (subtle) |
| **m7G** | 7-methylguanosine | G→A mismatches |
| **m5U** | 5-methyluridine | U mismatches (weak) |
| **Ψ** | pseudouridine | U→C mismatches (CMC-dependent) |
| **D** | dihydrouridine | U deletions/mismatches + RT stops |
| **I** | inosine | A→G mismatches |
| **i6A** | N6-isopentenyladenosine | A→G mismatches |
| **t6A** | N6-threonylcarbamoyladenosine | A mismatches |
| **m2G** | N2-methylguanosine | G→A mismatches |
| **m2₂G** | N2,N2-dimethylguanosine | G→A mismatches + RT stops |
| **acp3U** | 3-(3-amino-3-carboxypropyl)uridine | U mismatches |
| **Gm, Cm, Um, Am** | 2'-O-methylated nucleosides | Weak mismatches |
| **s4U** | 4-thiouridine | U→C mismatches |
| **mnm5s2U** | 5-methylaminomethyl-2-thiouridine | U mismatches + RT stops |
| **Q** | queuosine | G→any mismatches |
| **yW** | wybutosine | G mismatches + RT stops |
| **k2C** | lysidine | C→A mismatches |

## CLI Usage

### As part of the pipeline

```yaml
# config.yaml
run_modification_analysis: true
organism: "ecoli"           # ecoli, human, or mouse
no_modomics: false          # Set true to skip MODOMICS lookup
discover_novel_modifications: true
modification_min_coverage: 50
modification_alpha: 0.01
```

```bash
python -m trnaseq pipeline --config config.yaml --project-dir ./ --stages 6
```

### Standalone

```bash
python -m trnaseq modifications \
    --json-dir data/SWalign/ \
    --reference tRNA_database/ecoli/ecoli.fa \
    -o results/modifications/ \
    --organism ecoli \
    --discover-novel \
    --n-jobs 16
```

## Python API

```python
from trnaseq.modifications import (
    PositionalExtractor,
    RTSignatureAnalyzer,
    ModificationCaller,
    MODOMICSAnnotator,
)
from trnaseq.modifications.crosstalk import CrosstalkAnalyzer

# 1. Build PSCMs from alignment JSONs
extractor = PositionalExtractor(reference_fasta="ecoli.fa")
pscm_dict = extractor.run_parallel(json_paths, n_jobs=8)
# pscm_dict: {sample: {trna_name: numpy array (ref_len x 8)}}

# 2. Annotate known modifications from MODOMICS
annotator = MODOMICSAnnotator(organism="ecoli")
for trna_name, ref_seq in ref_sequences.items():
    known = annotator.get_known_mods_linear(trna_name, ref_seq)
    # known: [(position, short_name, full_name, rt_signature), ...]

# 3. Call modifications
caller = ModificationCaller(alpha=0.01)
# ... see modification_caller.py for full API

# 4. Crosstalk analysis
ct = CrosstalkAnalyzer(min_coverage=50, alpha=0.05)
ct_results = ct.analyze_multiple_samples(json_paths, mod_positions, n_jobs=8)
```

## Output Files

Stage 6 produces the following in `results/modifications/`:

| File | Description |
|------|-------------|
| `pscm_*.parquet` | Position-Specific Count Matrices per sample |
| `modification_calls.parquet` | Per-sample, per-position modification calls |
| `aggregated_calls.parquet` | Replicate-aggregated calls (majority vote) |
| `modification_summary.parquet` | Per-tRNA summary (positions, types, rates) |
| `crosstalk_results.parquet` | SLAC coordination analysis results |
| `modification_report.html` | Interactive dashboard |

## References

- Behrens et al. 2023 (NAR) — SLAC methodology for crosstalk analysis
- MODOMICS database — tRNA modification reference (https://genesilico.pl/modomics/)
- Wang et al. 2021 — RT stop methodology
