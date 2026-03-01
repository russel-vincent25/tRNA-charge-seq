# Charge tRNA-Seq

Pipeline for processing and analyzing charge tRNA-Seq data as described in our [manuscript](https://www.biorxiv.org/content/10.1101/2023.07.31.551363v1). Takes raw paired-end reads through adapter removal, barcode splitting, UMI deduplication, Smith-Waterman alignment, and produces tRNA charge quantification with interactive QC reports.

Only tested on Linux and macOS.

## Repository Structure

```
tRNA-charge-seq/
├── src/                        # Core preprocessing classes
│   ├── read_processing.py      #   Adapter removal, BC split, UMI trim
│   ├── alignment.py            #   SWIPE Smith-Waterman alignment
│   ├── stats_collection.py     #   Per-sample stats aggregation
│   ├── plotting.py             #   TRNA_plot visualization class
│   └── misc.py                 #   Utilities (read_tRNAdb_info, etc.)
│
├── trnaseq/                    # Analysis package + CLI
│   ├── cli/                    #   Command-line interface
│   │   ├── run_pipeline.py     #     Unified preprocessing pipeline
│   │   └── commands/           #     Subcommands (quantify, view, abundance)
│   ├── charge/                 #   Charge quantification (ChargeQuantifier)
│   ├── qc/                     #   QC dashboard (Plotly, 8 interactive panels)
│   ├── abundance/              #   Differential abundance (pyDESeq2)
│   ├── visualization/          #   Alignment viewer
│   ├── modifications/          #   RT signature & modification calling
│   ├── fragments/              #   Fragment classification & RT drop-off
│   └── io/                     #   Data storage (Parquet)
│
├── projects/                   # Analysis projects
│   ├── example/                #   Minimal example (8 samples, raw → plots)
│   └── RT-comp_script/         #   RT enzyme comparison (36 samples)
│
├── hpc/slurm/                  # SLURM job files for HMS O2
│   ├── submit_pipeline.sh      #   Master launcher (chains 3 jobs)
│   ├── stage0ab.job            #   Merge + barcode split
│   ├── stage0c_1.job           #   UMI + alignment (array, per sample)
│   └── stage2_6.job            #   Stats + charge + fragments + QC + mods
│
├── utils/                      # Index lists, scoring matrices
├── tRNA_database/              # Reference tRNA sequences
├── tRNA_database_masked/       # Masked references (for SWIPE)
├── tests/                      # Test suite
└── environment.yml             # Conda environment specification
```

## Installation

### 1. Create conda environment

```bash
conda env create -f environment.yml
conda activate tRNA-seq
```

This installs Python 3.10, all required packages (pandas, numpy, scipy, plotly, etc.), and bioinformatics tools (AdapterRemoval, SWIPE, BLAST).

### 2. Verify installation

```bash
python -c "import pandas, numpy, scipy, plotly; print('OK')"
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
```

> **Version constraints:** scipy <1.15 (seaborn compat), numpy <1.26 + pandas <2.2 (Excel reading compat)

## Two Workflows

### Notebook workflow (interactive analysis)

Best for exploring data, making plots, and small datasets.

1. Copy `projects/example/` to a new folder
2. Edit `sample_list.xlsx` with your samples
3. Run `process_data.ipynb` cell by cell

Classes are imported directly from `src/`:
```python
from src.read_processing import AR_merge, BC_split, UMI_trim
from src.alignment import SWIPE_align
from src.stats_collection import STATS_collection
from src.plotting import TRNA_plot
```

### CLI workflow (batch processing)

Best for large datasets and HPC deployment.

```bash
python -m trnaseq pipeline \
    --config config.yaml \
    --project-dir /path/to/project/ \
    --n-jobs 4
```

This runs the full pipeline (stages 0a–7). Progress is logged to `pipeline.log` in the project directory.

## Pipeline Stages

| Stage | What it does | Output |
|-------|-------------|--------|
| 0a | AdapterRemoval: merge paired reads | `data/AdapterRemoval/` |
| 0b | Barcode splitting | `data/BC_split/` |
| 0c | UMI trimming | `data/UMI_trimmed/` |
| 1 | SWIPE alignment | `data/SWalign/` |
| 2 | Stats collection + aggregation | `data/stats_collection/ALL_stats_aggregate.csv` |
| 3 | Charge quantification + fragment analysis | `charge_analysis/`, `fragment_analysis/` |
| 4 | Parquet storage (optional) | `parquet_data/` |
| 5 | QC dashboard | `qc_reports/QC_report.html`, `qc_reports/QC_summary.csv` |
| 6 | Modification analysis (optional) | `modification_analysis/` |
| 7 | Differential abundance (optional) | `abundance_analysis/` |

Run specific stages with `--stages 0a,0b` or `--stages 2,3,5`.

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
python -m trnaseq fragments --stats-dir data/stats_collection/ -o fragment_analysis/

# Modification analysis
python -m trnaseq modifications --json-dir data/SWalign/ \
    --reference tRNA_database/ecoli/ecoli.fa -o modification_analysis/ \
    --organism ecoli --discover-novel

# Differential abundance (requires pyDESeq2)
python -m trnaseq abundance -i ALL_stats_aggregate.csv \
    --sample-df sample_df.xlsx --level aa --control WT -o results/
```

## Configuration

Create a YAML config file (see `projects/example/config_beta_test.yaml` for a complete example):

```yaml
# Required
sample_list: "sample_list.xlsx"        # Relative to --project-dir
index_list: "utils/index_list.xlsx"    # Or absolute path
seq_dir: "raw_fastq"
tRNA_database:
  human: "/path/to/human-tRNAs.fa"
SWIPE_score_mat: "/path/to/nuc_score-matrix_2.txt"
common_seqs: null                      # Or path to common-seqs FASTA

# Processing
min_read_len: 39
min_score_align: 15
gap_penalty: 6
extension_penalty: 3
downsample_percentile: 50              # null for no downsampling
overwrite: true

# Stage 3: Charge + Fragments
run_charge_quantification: true
charge_count: "count"                  # "count" or "UMIcount"
charge_levels: [transcript, codon, aa]
include_mito: true
run_fragment_analysis: true

# Stage 4: Parquet (optional)
run_parquet_storage: false

# Stage 5: QC
run_qc_report: true

# Stage 6: Modifications (optional)
run_modification_analysis: false
organism: "ecoli"                      # ecoli, human, or mouse
discover_novel_modifications: false

# Stage 7: Abundance (optional, requires pyDESeq2)
run_abundance_analysis: false
# abundance_level: "aa"
# abundance_control: "WT"
```

Paths can be absolute or relative to `--project-dir`.

## HPC Deployment (SLURM)

For large datasets on HMS O2 or similar clusters. The pipeline supports per-sample parallelization via SLURM array jobs.

Set your project directory (must contain `data/raw_fastq/`, config, and sample list):

```bash
export PROJECT_DIR="/n/scratch/users/r/$USER/my_project"
```

Submit the pipeline from the **login node**:

```bash
# Automatic (chains 3 jobs with dependencies)
bash hpc/slurm/submit_pipeline.sh $PROJECT_DIR/config.yaml $PROJECT_DIR 72 32

# Or submit manually with sbatch
REPO="/home/$USER/github_repos/tRNA-charge-seq"

JOB0=$(sbatch --parsable $REPO/hpc/slurm/stage0ab.job $PROJECT_DIR/config.yaml $PROJECT_DIR)
JOB1=$(sbatch --parsable --dependency=afterok:$JOB0 --array=0-71%32 $REPO/hpc/slurm/stage0c_1.job $PROJECT_DIR/config.yaml $PROJECT_DIR)
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 $REPO/hpc/slurm/stage2_6.job $PROJECT_DIR/config.yaml $PROJECT_DIR)
```

See [`hpc/slurm/QUICK_START.md`](hpc/slurm/QUICK_START.md) for full setup and deployment guide.

## Input Data

- Raw paired-end reads in **bzip2-compressed FASTQ** format (`.fastq.bz2`)
- To convert from gzip: `ls *.gz | parallel "gunzip -c {} | bzip2 > {.}.bz2"`
- Sample list in Excel format with columns: `sample_name_unique`, `sample_name`, `replicate`, `fastq_mate1_filename`, `fastq_mate2_filename`, `P5_index`, `P7_index`, `barcode`, `species`

## QC Dashboard

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

## Dependencies

**Bioinformatics tools** (installed via conda):
- AdapterRemoval v2 — adapter trimming and read merging
- SWIPE — Smith-Waterman alignment
- BLAST (makeblastdb) — only for building new tRNA databases

**Key Python packages:**
- pandas, numpy, scipy, matplotlib, seaborn — data analysis
- plotly — interactive QC dashboard
- Biopython — sequence handling
- mpire — parallel processing
- pyDESeq2 — differential abundance (optional)

## Citation

If you use this pipeline, please cite:

> [Charge tRNA-Seq manuscript](https://www.biorxiv.org/content/10.1101/2023.07.31.551363v1)
