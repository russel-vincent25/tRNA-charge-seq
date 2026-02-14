# SLURM Pipeline Quick Start (HMS O2)

## Prerequisites

1. tRNA-charge-seq repo cloned to O2
2. Conda environment `tRNA-seq` with all dependencies
3. Config YAML pointing to correct paths on O2

## Usage

```bash
# From the repo root:
bash hpc/slurm/submit_pipeline.sh <config.yaml> <project_dir> [n_samples] [max_concurrent]
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `config.yaml` | Yes | Pipeline config (paths resolved relative to project_dir) |
| `project_dir` | Yes | Output directory for all pipeline results |
| `n_samples` | No | Number of samples (auto-detected from sample_list.xlsx) |
| `max_concurrent` | No | Max concurrent array tasks (default: 32) |

### Example (PJ39, 256 samples)

```bash
bash hpc/slurm/submit_pipeline.sh \
    /n/data/PJ39/config.yaml \
    /n/data/PJ39/ \
    256 \
    32
```

## Job Architecture

```
JOB0: stage0ab.sh         (single job)
  Stages 0a + 0b: Merge reads + BC split
  Resources: 4 CPU, 16G, 1 hour
      │
      ▼
JOB1: stage0c_1.sh        (array job, 0-N%32)
  Stages 0c + 1: UMI trim + SWIPE alignment (per sample)
  Resources: 4 CPU, 8G, 45 min per task
      │
      ▼
JOB2: stage2_5.sh         (single job)
  Stages 2 + 3 + 5: Stats + Charge + QC report
  Resources: 16 CPU, 64G, 2 hours
```

The `%32` throttle means at most 32 samples run simultaneously.
Adjust with the 4th argument if your admin prefers fewer concurrent jobs.

## Estimated Runtime (256 samples)

| Stage | Time | Notes |
|-------|------|-------|
| 0a+0b (JOB0) | ~35 min | Depends on number of input file pairs |
| 0c+1 (JOB1) | ~45 min | Wall time per task; all run in parallel |
| 2+3+5 (JOB2) | ~35 min | Sequential aggregation |
| **Total** | **~2 hours** | vs ~10 hours monolithic |

## Monitoring

```bash
# Check job status
squeue -u $USER

# Detailed job info
sacct -j <JOB0_ID>,<JOB1_ID>,<JOB2_ID> --format=JobID,State,Elapsed,MaxRSS

# Watch logs in real time
tail -f /path/to/PJ39/logs/*.out

# Check for failed array tasks
sacct -j <JOB1_ID> --format=JobID,State,ExitCode | grep -v COMPLETED
```

## Re-running Failed Tasks

If some array tasks fail (e.g., task 42 and 107):

```bash
# Re-run specific samples only
sbatch --array=42,107 hpc/slurm/stage0c_1.sh config.yaml /path/to/PJ39/

# Then re-run aggregation
sbatch hpc/slurm/stage2_5.sh config.yaml /path/to/PJ39/
```

## Resource Tuning

If jobs are getting killed (OOM) or timing out:

```bash
# Increase memory for alignment-heavy samples:
sbatch --mem=16G --array=0-255%32 hpc/slurm/stage0c_1.sh ...

# Increase time for large reference databases:
sbatch --time=01:30:00 --array=0-255%32 hpc/slurm/stage0c_1.sh ...

# Reduce concurrent tasks if admin requests:
sbatch --array=0-255%16 hpc/slurm/stage0c_1.sh ...
```

## Outputs

After completion, find results in `project_dir/`:

```
project_dir/
├── data/
│   ├── AdapterRemoval/    # Merged reads
│   ├── BC_split/          # Barcode-split FASTQs + read_length_distributions.csv
│   ├── UMI_trimmed/       # UMI-trimmed reads
│   ├── SWalign/           # Alignment JSONs
│   └── stats_collection/  # Per-sample + aggregate stats CSVs
├── charge_analysis/       # charge_df_{aa,codon,transcript}.csv
├── QC_summary.csv         # One row per sample, all metrics
├── QC_report.html         # Interactive Plotly dashboard
├── sample_df.xlsx         # Updated sample metadata
├── logs/                  # SLURM job logs
└── preprocessing.log      # Pipeline log
```
