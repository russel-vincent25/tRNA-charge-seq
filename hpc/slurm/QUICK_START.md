# SLURM Pipeline Deployment Guide (HMS O2)

## 1. Initial Setup

### Clone the repo

```bash
mkdir -p /home/$USER/github_repos
cd /home/$USER/github_repos

git clone https://github.com/russel-vincent25/tRNA-charge-seq.git
cd tRNA-charge-seq
```

### Create conda environment

Run on an interactive compute node (not the login node):

```bash
srun -pty -p interactive --mem 2G -t 0-01:00:00 /bin/bash

module load conda/miniforge3/24.11.3-0
conda env create -f environment.yml    # ~10-15 min
conda activate tRNA-seq

# Install trnaseq package (editable — survives git pulls)
pip install -e .

# Verify
python -c "import pandas, numpy, scipy, plotly, trnaseq; print('OK')"
AdapterRemoval --version
swipe -h 2>&1 | head -1
```

If the conda solver has issues, create manually:

```bash
conda create -n tRNA-seq python=3.10 -y
conda activate tRNA-seq
conda install -c conda-forge -c bioconda numpy">=1.24,<1.26" pandas">=2.0,<2.2" scipy">=1.11,<1.15" matplotlib seaborn==0.13.0 plotly logomaker biopython openpyxl xlrd pyarrow pyyaml mpire jellyfish json_stream natsort tqdm adapterremoval swipe blast bzip2 pigz pillow wand imagemagick pytest -y
pip install pydeseq2
pip install -e .
```

---

## 2. Project Directory Setup

The pipeline uses `--project-dir` as the working directory. It expects:

```
project_dir/
├── config.yaml              # You provide: pipeline config
├── sample_list.xlsx         # You provide: sample metadata
└── data/
    └── raw_fastq/           # You provide: *.fastq.bz2 files (R1 + R2)
```

The pipeline creates all output alongside `data/`:

```
project_dir/
├── config.yaml
├── sample_list.xlsx
├── sample_df.xlsx                   # Updated sample metadata
├── logs/
│   ├── pipeline.log                 # Pipeline log
│   └── computing_metrics.csv        # Per-stage timing
├── data/
│   ├── raw_fastq/                   # Original FASTQs (already here)
│   ├── AdapterRemoval/              # Merged reads
│   ├── BC_split/                    # Barcode-split FASTQs
│   ├── UMI_trimmed/                 # UMI-trimmed reads
│   ├── SWalign/                     # Alignment JSONs
│   └── stats_collection/            # Per-sample + aggregate stats
├── results/
│   ├── charge/                      # charge_df_{aa,codon,transcript}.csv
│   ├── fragments/                   # Fragment counts, RT drop-off, lengths
│   ├── parquet/                     # Parquet-format data (optional)
│   ├── modifications/               # PSCM, modification calls, crosstalk
│   └── abundance/                   # DESeq2 results + volcano/MA plots
└── qc_reports/                      # QC_report.html + QC_summary.csv
```

### Path resolution

Config paths can be **relative** (resolved against `--project-dir`) or **absolute**:

```yaml
# Relative — resolved to {project_dir}/sample_list.xlsx
sample_list: "sample_list.xlsx"

# Absolute — used as-is (recommended for repo-level files)
index_list: "/home/ruv988/github_repos/tRNA-charge-seq/utils/index_list_updated.xlsx"
```

Fields that support relative paths: `sample_list`, `index_list`, `SWIPE_score_mat`,
`common_seqs`, `adapter_sequences`, and all `tRNA_database` entries.

### If your data already has the right structure

Point `--project-dir` at the existing experiment directory and copy configs into it:

```bash
export PROJECT_DIR="/n/scratch/users/r/$USER/experiment/2025-01-24-run"
cp config.yaml sample_list.xlsx "$PROJECT_DIR/"
```

### If starting fresh

```bash
export PROJECT_DIR="/n/scratch/users/r/$USER/my_project"
mkdir -p "$PROJECT_DIR/data/raw_fastq"
# Copy or symlink FASTQs (R1 + R2 only, not I1 index reads)
ln -s /path/to/original/*_R1.fastq.bz2 "$PROJECT_DIR/data/raw_fastq/"
ln -s /path/to/original/*_R2.fastq.bz2 "$PROJECT_DIR/data/raw_fastq/"
cp config.yaml sample_list.xlsx "$PROJECT_DIR/"
```

### Verify before running

```bash
ls "$PROJECT_DIR/data/raw_fastq/"*_R1.fastq.bz2 | wc -l   # Should match N_SAMPLES
ls "$PROJECT_DIR/config"*.yaml "$PROJECT_DIR/sample_list"*  # Configs present
mkdir -p /home/$USER/jobOutput                               # For SLURM logs

# Run preflight checks
cd /home/$USER/github_repos/tRNA-charge-seq
python -m trnaseq pipeline --config "$PROJECT_DIR/config.yaml" --project-dir "$PROJECT_DIR" --preflight
```

---

## 3. How the Pipeline Works

Three SLURM `.job` files are chained with dependencies:

```
JOB0: stage0ab.job        (single job)
  Stages 0a + 0b: Merge reads + barcode split
      │
      ▼
JOB1: stage0c_1.job       (array job, 0-N samples)
  Stages 0c + 1: UMI trim + SWIPE alignment (per sample)
      │
      ▼
JOB2: stage2_6.job        (single job)
  Stages 2 + 3 + 5 + 6 + 7: Stats + Charge + Fragments + QC + Mods + Abundance
```

| File | Type | Default Resources |
|------|------|-------------------|
| `stage0ab.job` | Single job | Auto-scaled CPUs, 4G |
| `stage0c_1.job` | Array job | 2 CPU, 2G per task, 8h |
| `stage2_6.job` | Single job | 16 CPU, 16-32G, 2-12h |

**Thread awareness:** The pipeline uses `--threads-per-job` to control how many
threads each AdapterRemoval or SWIPE subprocess uses. Stage 0ab computes
`n_jobs = CPUs / threads_per_job`; stage 0c+1 gives all CPUs to a single SWIPE
process. The launcher auto-sizes all resources based on sample count.

Job logs go to `/home/$USER/jobOutput/tseq-*.out` and `*.err`.
Email notifications (BEGIN, END, FAIL) are sent for each job.

---

## 4. Running the Pipeline

### Option A: Automatic submission

The launcher script submits all 3 jobs from the **login node**:

```bash
cd /home/$USER/github_repos/tRNA-charge-seq

bash hpc/slurm/submit_pipeline.sh <config.yaml> <project_dir> [n_samples] [max_concurrent]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `config.yaml` | Yes | Path to pipeline config |
| `project_dir` | Yes | Path to project directory |
| `n_samples` | No | Number of samples (auto-detected if omitted) |
| `max_concurrent` | No | Max concurrent array tasks (default: 32) |

Example:

```bash
bash hpc/slurm/submit_pipeline.sh "$PROJECT_DIR/config.yaml" "$PROJECT_DIR"
```

### Option B: Manual submission

Submit each `.job` file yourself with `sbatch`:

```bash
REPO="/home/$USER/github_repos/tRNA-charge-seq"
CONFIG="$PROJECT_DIR/config.yaml"
N=71  # Number of samples minus 1 (0-indexed)

# Job 0: Merge + BC split
JOB0=$(sbatch --parsable \
    --cpus-per-task=16 --mem=4G -t 02:00:00 \
    "$REPO/hpc/slurm/stage0ab.job" "$CONFIG" "$PROJECT_DIR" 2)
echo "JOB0: $JOB0"

# Job 1: Per-sample UMI + alignment
JOB1=$(sbatch --parsable \
    --dependency=afterok:${JOB0} --array=0-${N}%32 \
    --cpus-per-task=2 --mem=2G -t 08:00:00 \
    "$REPO/hpc/slurm/stage0c_1.job" "$CONFIG" "$PROJECT_DIR")
echo "JOB1: $JOB1"

# Job 2: Stats + charge + fragments + QC + mods + abundance
JOB2=$(sbatch --parsable \
    --dependency=afterok:${JOB1} \
    --cpus-per-task=16 --mem=32G -t 04:00:00 \
    "$REPO/hpc/slurm/stage2_6.job" "$CONFIG" "$PROJECT_DIR")
echo "JOB2: $JOB2"
```

---

## 5. Estimated Runtimes

| Dataset | 0a+0b | 0c+1 (wall) | 2+3+5+6+7 | Total |
|---------|-------|-------------|------------|-------|
| 24 samples | ~15 min | ~59 min | ~50 min | ~2 hours |
| 72 samples | ~25 min | ~59 min | ~2 hr | ~3 hours |
| 264 samples | ~35 min | ~59 min | ~5 hr | ~6 hours |

Stage 2+3+5+6+7 includes ~30 min overhead for modification analysis (PSCM extraction + SLAC crosstalk). The launcher allocates a minimum of 2 hours.

---

## 6. Monitoring

```bash
# Check job queue
squeue -u $USER

# Detailed job info (replace with your job IDs)
sacct -j <JOB0>,<JOB1>,<JOB2> --format=JobID,JobName,State,Elapsed,MaxRSS

# Watch logs
tail -f /home/$USER/jobOutput/tseq-*.out

# Check for failed array tasks
sacct -j <JOB1> --format=JobID,State,ExitCode | grep -v COMPLETED
```

---

## 7. Verifying Results

```bash
# Check output directories
ls "$PROJECT_DIR/qc_reports/QC_report.html" \
   "$PROJECT_DIR/results/charge/" \
   "$PROJECT_DIR/results/fragments/" \
   "$PROJECT_DIR/results/modifications/"

# Quick sanity check
python -c "
import pandas as pd
qc = pd.read_csv('$PROJECT_DIR/qc_reports/QC_summary.csv')
print(qc[['sample_name_unique', 'N_pairs', 'N_total_aligned']].to_string())
print(f'\n{len(qc)} samples processed')
"

# Download QC report to view in browser
scp $USER@transfer.rc.hms.harvard.edu:$PROJECT_DIR/qc_reports/QC_report.html ~/Desktop/
```

---

## 8. Re-running Failed Tasks

If specific array tasks fail (e.g., tasks 5 and 12):

```bash
# Re-run only failed samples
sbatch --array=5,12 $REPO/hpc/slurm/stage0c_1.job "$CONFIG" "$PROJECT_DIR"

# Then re-run aggregation
sbatch $REPO/hpc/slurm/stage2_6.job "$CONFIG" "$PROJECT_DIR"
```

To re-run only specific analysis stages (e.g., just modifications):

```bash
sbatch --mem=32G -c 16 -p short -t 2:00:00 \
    --wrap "module load conda/miniforge3/24.11.3-0; conda activate tRNA-seq; \
    python -m trnaseq pipeline --config $CONFIG --project-dir $PROJECT_DIR --stages 6"
```

---

## 9. Resource Tuning

The launcher (`submit_pipeline.sh`) auto-computes CPUs, memory, and wall time
from `N_SAMPLES`. For manual runs, override `.job` defaults on the command line:

```bash
# More memory for large sample counts
sbatch --mem=64G $REPO/hpc/slurm/stage2_6.job "$CONFIG" "$PROJECT_DIR"

# More time for modification analysis with many samples
sbatch -t 08:00:00 $REPO/hpc/slurm/stage2_6.job "$CONFIG" "$PROJECT_DIR"

# Fewer concurrent tasks (if admin requests)
sbatch --array=0-71%16 $REPO/hpc/slurm/stage0c_1.job "$CONFIG" "$PROJECT_DIR"
```

**Thread tuning:** `threads_per_job` controls how many threads each AdapterRemoval
or SWIPE subprocess uses. The pipeline computes `n_jobs = CPUs / threads_per_job`.
Set via config YAML (`threads_per_job: 2`) or CLI (`--threads-per-job 2`).
Default is 2.

---

## 10. Clean Re-run

To start fresh (e.g., after config changes or code updates):

```bash
cd $PROJECT_DIR

# Remove all outputs
rm -rf results/ qc_reports/ logs/
rm -f inp_file_df.xlsx sample_df.xlsx

# Remove per-sample data (forces full reprocess)
rm -rf data/AdapterRemoval data/BC_split data/UMI_trimmed data/SWalign data/stats_collection

# Pull latest code and reinstall
cd /home/$USER/github_repos/tRNA-charge-seq
git pull
pip install -e .

# Resubmit
bash hpc/slurm/submit_pipeline.sh "$PROJECT_DIR/config.yaml" "$PROJECT_DIR"
```

---

## 11. Troubleshooting

| Problem | Solution |
|---------|----------|
| `conda: command not found` | `module load conda/miniforge3/24.11.3-0` |
| `No module named trnaseq` | `pip install -e .` from the repo root |
| `requires a different Python: 3.9` | Supported — `requires-python >= 3.9` |
| `ModuleNotFoundError: pydeseq2` | `pip install pydeseq2` |
| `GLIBCXX_3.4.30 not found` (scipy) | `pip install scipy==1.11.4` |
| Array task OOM killed | Resubmit with `--mem=16G` |
| Array task timeout | Resubmit with `-t 01:30:00` |
| Stage 2 OOM with many samples | Use `--mem=64G` or reduce `--n-jobs 4` |
| `FileNotFoundError: sample_list` | Check sample list is in `$PROJECT_DIR` |
| `FileNotFoundError: raw_fastq` | Check path: `ls $PROJECT_DIR/data/raw_fastq/` |
| Stage 2 fails (no stats) | Check JOB1 logs — some array tasks may have failed |
| `tRNA_database` not found | Verify absolute path in config exists on O2 |
| `swipe: command not found` | Conda env not activated — check `.job` module load |
