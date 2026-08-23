# tRNA-charge-seq

Pipeline for processing and analyzing charge tRNA-Seq data as described in our [manuscript](https://www.biorxiv.org/content/10.1101/2023.07.31.551363v1). Takes raw paired-end reads through adapter removal, barcode splitting, UMI deduplication, Smith-Waterman alignment, and produces tRNA charge quantification, modification calling, and interactive QC reports.

Only tested on Linux and macOS.

## Repository Structure

```
tRNA-charge-seq/
├── trnaseq/                    # Analysis package + CLI
│   ├── cli/                    #   Command-line interface
│   │   ├── run_pipeline.py     #     Unified preprocessing pipeline
│   │   └── commands/           #     Subcommands (quantify, view, abundance, ...)
│   ├── charge/                 #   Charge quantification (ChargeQuantifier)
│   ├── qc/                     #   QC + report dashboards (Plotly, interactive)
│   ├── modifications/          #   Modification analysis (PSCM, MODOMICS, crosstalk)
│   ├── fragments/              #   Fragment classification & RT drop-off
│   ├── abundance/              #   Differential abundance (pyDESeq2)
│   ├── visualization/          #   Alignment viewer
│   └── io/                     #   Data storage (Parquet)
│
├── src/                        # Core preprocessing classes
│   ├── read_processing.py      #   Adapter removal, BC split, UMI trim
│   ├── alignment.py            #   SWIPE Smith-Waterman alignment
│   ├── stats_collection.py     #   Per-sample stats aggregation
│   ├── plotting.py             #   TRNA_plot visualization class
│   └── misc.py                 #   Utilities (read_tRNAdb_info, etc.)
│
├── hpc/slurm/                  # SLURM job files for HMS O2
│   ├── submit_pipeline.sh      #   Master launcher (chains 3 jobs)
│   ├── stage0ab.job            #   Merge + barcode split
│   ├── stage0c_1.job           #   UMI + alignment (array, per sample)
│   └── stage2_6.job            #   Stats + charge + fragments + QC + mods + abundance
│
├── projects/                   # Analysis projects
│   └── example_script_test/    #   Example config + sample list
│
├── utils/                      # Index lists, scoring matrices, adapter sequences
├── tRNA_database/              # Reference tRNA sequences (FASTA + BLAST dbs)
├── tests/                      # Test suite
├── pyproject.toml              # Package metadata
└── environment.yml             # Conda environment specification
```

## Installation

### 1. Create conda environment

```bash
conda env create -f environment.yml
conda activate tRNA-seq
```

### 2. Install the trnaseq package

```bash
pip install -e .
```

This installs `trnaseq` in editable mode — survives `git pull` without reinstalling.

### 3. Verify installation

```bash
python -c "import trnaseq; print('OK')"
AdapterRemoval --version
swipe -h 2>&1 | head -1
```

### Manual install (if conda solver has issues)

```bash
conda create -n tRNA-seq python=3.10 -y
conda activate tRNA-seq
conda install -c conda-forge -c bioconda \
    numpy">=1.24,<1.26" pandas">=2.0,<2.2" scipy">=1.11,<1.15" \
    matplotlib seaborn==0.13.0 plotly logomaker \
    biopython openpyxl xlrd pyarrow pyyaml \
    mpire jellyfish json_stream natsort tqdm \
    adapterremoval swipe blast bzip2 pigz \
    pillow wand imagemagick pytest -y
pip install pydeseq2
pip install -e .
```

> **Version constraints:** scipy <1.15 (seaborn compat), numpy <1.26 + pandas <2.2 (Excel reading compat). Python >=3.9 supported.

## Quick Start

```bash
python -m trnaseq pipeline \
    --config config.yaml \
    --project-dir /path/to/project/ \
    --n-jobs 4
```

This runs the full pipeline (stages 0a–7). Progress is logged to `logs/pipeline.log` in the project directory.

### Preflight validation

Before submitting long-running jobs, verify your setup:

```bash
python -m trnaseq pipeline --config config.yaml --project-dir ./ --preflight
```

This checks that all input files, databases, and tools are accessible without processing any data.

## Pipeline Stages

| Stage | What it does | Output |
|-------|-------------|--------|
| 0a | AdapterRemoval: merge paired reads | `data/AdapterRemoval/` |
| 0b | Barcode splitting | `data/BC_split/` |
| 0c | UMI trimming (anchored or pyrimidine mode) | `data/UMI_trimmed/` |
| 1 | SWIPE Smith-Waterman alignment | `data/SWalign/` |
| 2 | Stats collection + aggregation | `data/stats_collection/` |
| 3 | Charge quantification + fragment analysis | `results/charge/`, `results/fragments/` |
| 4 | Parquet storage (optional) | `results/parquet/` |
| 5 | QC dashboard | `qc_reports/` |
| 6 | Modification analysis (optional) | `results/modifications/` |
| 7 | Differential abundance (optional) | `results/abundance/` |

Run specific stages: `--stages 0a,0b` or `--stages 2,3,5`.

## Configuration

Create a YAML config file (see `projects/example_script_test/config.yaml` for a complete example):

```yaml
# Required
sample_list: "sample_list.xlsx"        # Relative to --project-dir
index_list: "/path/to/index_list.xlsx" # Or absolute path
seq_dir: "raw_fastq"
tRNA_database:
  human: "/path/to/human-tRNAs.fa"
SWIPE_score_mat: "/path/to/nuc_score-matrix_2.txt"
common_seqs: null                      # Or path to common-seqs FASTA.
                                       # MUST be null for --sample-index (array) runs: it
                                       # decompresses to a fixed path in the pipeline repo, so
                                       # concurrent tasks race on one file, and _collect_stats
                                       # then overwrites the real N_mapped with the common-seq
                                       # count. Only set it if the FASTA actually matches your
                                       # library -- check *_common-seq-obs.json is not all zeros.
realign_overwrite: true                # false = stage 1 reuses an existing *_SWalign.json.bz2
                                       # instead of realigning. Does NOT check the reference,
                                       # score matrix, or reads still match -- opt in per run.

# Read processing
min_read_len: 39
min_score_align: 15
gap_penalty: 6
extension_penalty: 3
overwrite: true
threads_per_job: 2

# UMI trimming
umi_trim_mode: "anchored"              # "anchored" (GCv4) or "pyrimidine" (legacy)
adapter_sequences: "/path/to/adapter_sequences.yaml"
umi_anchor: "GCv4"                     # Anchor name (anchored mode only)
umi_max_stagger: 3
umi_anchor_max_dist: 1

# Downsampling
downsample_percentile: 50              # null for no downsampling
downsample_absolute: null              # Or fixed read count

# Stage 3: Charge + Fragments
run_charge_quantification: true
charge_count: "count"                  # "count" or "UMIcount"
charge_levels: [transcript, codon, aa]
include_synthetic: false
include_mito: true
run_fragment_analysis: true
fragment_min_reads: 10
fragment_write_csv: true               # Also write CSV (default: Parquet only)

# Stage 4: Parquet (optional)
run_parquet_storage: false

# Stage 5: QC
run_qc_report: true

# Stage 6: Modifications (optional)
run_modification_analysis: false
organism: "ecoli"                      # ecoli, human, or mouse
no_modomics: false                     # true to skip MODOMICS lookup
discover_novel_modifications: false
modification_min_coverage: 50
modification_alpha: 0.01

# Stage 7: Abundance (optional)
run_abundance_analysis: false
abundance_level: "aa"                  # transcript, codon, or aa
abundance_control: "WT"               # Control condition name
```

Paths can be absolute or relative to `--project-dir`.

## CLI Commands

```bash
# Full pipeline
python -m trnaseq pipeline --config config.yaml --project-dir ./

# Charge quantification only
python -m trnaseq quantify -i ALL_stats_aggregate.csv -o charge.csv

# Alignment viewer
python -m trnaseq view --json data/SWalign/sample.json.bz2 --trna tRNA-Ala-AGC-1-1
python -m trnaseq view --json data/SWalign/sample.json.bz2 --list

# Fragment analysis
python -m trnaseq fragments --stats-dir data/stats_collection/ -o results/fragments/

# Modification analysis
python -m trnaseq modifications --json-dir data/SWalign/ \
    --reference tRNA_database/ecoli/ecoli.fa -o results/modifications/ \
    --organism ecoli --discover-novel

# Differential abundance (requires pyDESeq2)
python -m trnaseq abundance -i ALL_stats_aggregate.csv \
    --sample-df sample_df.xlsx --level aa --control WT -o results/

# Build MODOMICS reference FASTA
python -m trnaseq build-reference --organism ecoli -o ecoli_modomics_ref.fa
```

## Output Directory Layout

After a full pipeline run, the project directory contains:

```
project_dir/
├── config.yaml                      # Your config
├── sample_list.xlsx                 # Your sample metadata
├── sample_df.xlsx                   # Updated sample metadata (pipeline-generated)
├── logs/
│   ├── pipeline.log                 # Pipeline log
│   └── computing_metrics.csv        # Per-stage timing
├── data/
│   ├── raw_fastq/                   # Input FASTQs
│   ├── AdapterRemoval/              # Merged reads
│   ├── BC_split/                    # Barcode-split FASTQs
│   ├── UMI_trimmed/                 # UMI-trimmed reads
│   ├── SWalign/                     # Alignment JSONs (*.json.bz2)
│   └── stats_collection/            # Per-sample stats + ALL_stats_aggregate.csv
├── results/
│   ├── charge/                      # charge_df_{aa,codon,transcript}.csv
│   ├── fragments/                   # Fragment counts, RT drop-off, lengths
│   ├── parquet/                     # Parquet-format data (optional)
│   ├── modifications/               # Modification calls, PSCM, crosstalk
│   └── abundance/                   # DESeq2 results + reports
└── qc_reports/                      # QC_report.html + QC_summary.csv
```

## Modification Analysis (Stage 6, optional)

Stage 6 runs a 7-phase modification analysis pipeline. Disabled by default — enable with `run_modification_analysis: true`.

1. **PSCM extraction** — Builds Position-Specific Count Matrices from SWalign JSONs (parallelized)
2. **Background estimation** — Estimates per-sample error rates from low-modification positions
3. **Modification calling** — Tests each position against background using binomial test + BH-FDR
4. **Replicate aggregation** — Merges calls across replicates (majority vote + meta-analysis)
5. **Summary** — Generates per-tRNA modification summary tables
6. **SLAC crosstalk** — Analyzes per-read modification coordination (Fisher's exact test, parallelized)
7. **Report** — Interactive HTML dashboard

MODOMICS database integration provides annotation of known modifications with alignment-based position mapping. When MODOMICS data is unavailable for a specific tRNA, isotype-level modification transfer is used (excluding anticodon loop positions 32–38).

23 modification profiles are supported, including m1A, m1G, m3C, m5C, m7G, Ψ, D, I, i6A, t6A, m2G, m2₂G, and more.

See [`trnaseq/modifications/README.md`](trnaseq/modifications/README.md) for API details.

## Differential Abundance Analysis (Stage 7, optional)

Stage 7 performs differential tRNA abundance analysis using [pyDESeq2](https://pydeseq2.readthedocs.io/). It compares tRNA expression levels between experimental conditions and a control. Disabled by default — enable with `run_abundance_analysis: true`.

```yaml
# config.yaml
run_abundance_analysis: true
abundance_level: "aa"          # "transcript", "codon", or "aa"
abundance_control: "WT"        # Control condition (must match sample_name in sample list)
```

- **Levels**: Analyze at transcript (individual tRNA), codon (anticodon), or amino acid level
- **Statistical method**: DESeq2 negative binomial model with Wald test
- **Outputs** (`results/abundance/`):
  - DESeq2 results table (log2 fold change, p-value, adjusted p-value)
  - Interactive HTML report with volcano plot, MA plot, DE summary table, and PCA

Standalone usage:

```bash
python -m trnaseq abundance \
    -i data/stats_collection/ALL_stats_aggregate.csv \
    --sample-df sample_df.xlsx \
    --level aa \
    --control WT \
    -o results/abundance/
```

## QC Dashboard (Stage 5)

Stage 5 generates an interactive HTML report (`qc_reports/QC_report.html`) with 8 panels:

1. **Summary table** — per-sample metrics with pass/warn/fail highlighting
2. **Pipeline funnel** — read counts through each processing stage
3. **Mapping rate** — alignment success per sample
4. **UMI saturation** — unique vs total UMI counts
5. **Fragment types** — stacked bar (full-length, RT drop-off, 5'tRF, degraded)
6. **Read length distribution** — overlapped density per sample
7. **PCA** — replicate clustering colored by sample group
8. **Replicate correlation** — log-scale RPM scatter with Pearson R

Self-contained (no internet required), uses Plotly for interactivity.

## HPC Deployment (SLURM)

For large datasets on HMS O2 or similar clusters:

```bash
bash hpc/slurm/submit_pipeline.sh config.yaml /path/to/project_dir
```

This auto-detects sample count and submits 3 chained jobs:
- **JOB0** (stage0ab): Merge + barcode split
- **JOB1** (stage0c_1, array): UMI + alignment (1 task per sample)
- **JOB2** (stage2_6): Stats + charge + fragments + QC + modifications + abundance

Resources are auto-scaled based on sample count. See [`hpc/slurm/QUICK_START.md`](hpc/slurm/QUICK_START.md) for the full deployment guide.

## Notebook Workflow (interactive)

For exploring data interactively with small datasets:

1. Copy `projects/example_script_test/` to a new folder
2. Edit `sample_list.xlsx` with your samples
3. Import classes directly from `src/`:

```python
from src.read_processing import AR_merge, BC_split, UMI_trim
from src.alignment import SWIPE_align
from src.stats_collection import STATS_collection
from src.plotting import TRNA_plot
```

## Input Data

- Raw paired-end reads in **bzip2-compressed FASTQ** format (`.fastq.bz2`)
- To convert from gzip: `ls *.gz | parallel "gunzip -c {} | bzip2 > {.}.bz2"`
- Sample list in Excel format with columns: `sample_name_unique`, `sample_name`, `replicate`, `fastq_mate1_filename`, `fastq_mate2_filename`, `P5_index`, `P7_index`, `barcode`, `species`

## Dependencies

**Bioinformatics tools** (installed via conda):
- AdapterRemoval v2 — adapter trimming and read merging
- SWIPE — Smith-Waterman alignment
- BLAST (makeblastdb) — only for building new tRNA databases

**Key Python packages:**
- pandas, numpy, scipy, matplotlib, seaborn — data analysis
- plotly — interactive QC dashboard
- Biopython — sequence handling and pairwise alignment
- mpire — parallel processing
- pyDESeq2 — differential abundance analysis

## Citation

If you use this pipeline, please cite:

> [Charge tRNA-Seq manuscript](https://www.biorxiv.org/content/10.1101/2023.07.31.551363v1)
