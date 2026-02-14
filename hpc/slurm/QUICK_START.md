# SLURM Pipeline Quick Start (HMS O2)

## Prerequisites

1. tRNA-charge-seq repo cloned to O2
2. Conda environment `tRNA-seq` with all dependencies (see `environment.yml`)
3. Config YAML pointing to correct paths on O2

## How It Works

The pipeline uses **3 SLURM job files** chained with dependencies.
You can either use the `submit_pipeline.sh` launcher (submits all 3 automatically)
or submit each `.job` file manually with `sbatch`.

### Job Files

| File | Type | What it does | Resources |
|------|------|-------------|-----------|
| `stage0ab.job` | Single job | Merge reads + barcode split | 4 CPU, 16G, 1 hr |
| `stage0c_1.job` | Array job | UMI trim + SWIPE alignment (per sample) | 4 CPU, 8G, 45 min |
| `stage2_5.job` | Single job | Stats + charge + QC report | 16 CPU, 64G, 2 hr |

All jobs log to `/home/ruv988/jobOutput/` and send email notifications
(BEGIN, END, FAIL) to `russel_vincent@hms.harvard.edu`.

### Job Chain

```
JOB0: stage0ab.job        (single job)
  Stages 0a + 0b: Merge reads + BC split
      │
      ▼
JOB1: stage0c_1.job       (array job, 0-N%32)
  Stages 0c + 1: UMI trim + SWIPE alignment (per sample)
      │
      ▼
JOB2: stage2_5.job        (single job)
  Stages 2 + 3 + 5: Stats + Charge + QC report
```

## Option A: Automatic Submission

```bash
# From the repo root on the login node:
bash hpc/slurm/submit_pipeline.sh <config.yaml> <project_dir> [n_samples] [max_concurrent]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `config.yaml` | Yes | Pipeline config (paths resolved relative to project_dir) |
| `project_dir` | Yes | Output directory for all pipeline results |
| `n_samples` | No | Number of samples (auto-detected from sample_list.xlsx) |
| `max_concurrent` | No | Max concurrent array tasks (default: 32) |

Example (YLDC, 72 samples):
```bash
bash hpc/slurm/submit_pipeline.sh \
    /home/ruv988/projects/YLDC/config_YLDC.yaml \
    /home/ruv988/projects/YLDC/ \
    72 32
```

## Option B: Manual Submission

Submit each job yourself with `sbatch`:

```bash
REPO="/home/ruv988/github_repos/tRNA-charge-seq"
CONFIG="/home/ruv988/projects/YLDC/config_YLDC.yaml"
PROJECT_DIR="/home/ruv988/projects/YLDC"

# Step 1: Merge + BC split
JOB0=$(sbatch --parsable "$REPO/hpc/slurm/stage0ab.job" "$CONFIG" "$PROJECT_DIR")
echo "JOB0: $JOB0"

# Step 2: Per-sample alignment (wait for JOB0)
JOB1=$(sbatch --parsable \
    --dependency=afterok:${JOB0} \
    --array=0-71%32 \
    "$REPO/hpc/slurm/stage0c_1.job" "$CONFIG" "$PROJECT_DIR")
echo "JOB1: $JOB1"

# Step 3: Aggregation (wait for all JOB1 tasks)
JOB2=$(sbatch --parsable \
    --dependency=afterok:${JOB1} \
    "$REPO/hpc/slurm/stage2_5.job" "$CONFIG" "$PROJECT_DIR")
echo "JOB2: $JOB2"
```

## Estimated Runtimes

| Dataset | 0a+0b | 0c+1 (wall) | 2+3+5 | Total |
|---------|-------|-------------|-------|-------|
| 24 samples | ~15 min | ~5 min | ~10 min | **~30 min** |
| 72 samples | ~25 min | ~10 min | ~15 min | **~50 min** |
| 256 samples | ~35 min | ~45 min | ~35 min | **~2 hours** |

## Monitoring

```bash
# Check job queue
squeue -u $USER

# Detailed job info
sacct -j <JOB0>,<JOB1>,<JOB2> --format=JobID,JobName,State,Elapsed,MaxRSS

# Watch logs
tail -f /home/ruv988/jobOutput/tseq-*.out

# Check for failed array tasks
sacct -j <JOB1> --format=JobID,State,ExitCode | grep -v COMPLETED
```

## Re-running Failed Tasks

If some array tasks fail (e.g., tasks 5 and 12):

```bash
# Re-run specific samples only
sbatch --array=5,12 hpc/slurm/stage0c_1.job $CONFIG $PROJECT_DIR

# Then re-run aggregation
sbatch hpc/slurm/stage2_5.job $CONFIG $PROJECT_DIR
```

## Resource Tuning

If jobs are getting OOM-killed or timing out, override defaults on the command line:

```bash
# More memory for alignment-heavy samples
sbatch --mem=16G --array=0-71%32 hpc/slurm/stage0c_1.job ...

# More time
sbatch -t 01:30:00 --array=0-71%32 hpc/slurm/stage0c_1.job ...

# Fewer concurrent tasks
sbatch --array=0-71%16 hpc/slurm/stage0c_1.job ...
```

## Outputs

After completion, find results in `project_dir/`:

```
project_dir/
├── ALL_stats_aggregate.csv   # All tRNA counts per sample
├── QC_report.html            # Interactive Plotly dashboard
├── QC_summary.csv            # One row per sample, all metrics
├── charge_analysis/          # charge_df_{aa,codon,transcript}.csv
├── sample_df.xlsx            # Updated sample metadata
├── preprocessing.log         # Pipeline log
└── data/
    ├── AdapterRemoval/       # Merged reads
    ├── BC_split/             # Barcode-split FASTQs
    ├── UMI_trimmed/          # UMI-trimmed reads
    ├── SWalign/              # Alignment JSONs
    └── stats_collection/     # Per-sample + aggregate stats
```
