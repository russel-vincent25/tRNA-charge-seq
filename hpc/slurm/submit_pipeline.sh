#!/bin/bash
# ==============================================================================
# Master launcher for tRNA-charge-seq pipeline on HMS O2 (SLURM)
# ==============================================================================
#
# This script runs on the LOGIN NODE. It submits 3 sbatch jobs that execute
# on compute nodes with dependency chaining:
#
#   JOB0 (stage0ab.job)  →  JOB1 (stage0c_1.job, array)  →  JOB2 (stage2_5.job)
#
# Usage:
#   bash submit_pipeline.sh <config.yaml> <project_dir> [n_samples] [max_concurrent]
#
# Arguments:
#   config.yaml     - Pipeline configuration file
#   project_dir     - Project output directory
#   n_samples       - Number of samples (auto-detected from sample_list if omitted)
#   max_concurrent  - Max concurrent array tasks (default: 32)
#
# Example:
#   bash submit_pipeline.sh /home/ruv988/projects/YLDC/config_YLDC.yaml \
#                           /home/ruv988/projects/YLDC/ 72 32
#
# ==============================================================================

set -euo pipefail

CONFIG="${1:?Usage: submit_pipeline.sh <config.yaml> <project_dir> [n_samples] [max_concurrent]}"
PROJECT_DIR="${2:?Usage: submit_pipeline.sh <config.yaml> <project_dir> [n_samples] [max_concurrent]}"
MAX_CONCURRENT="${4:-32}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Resolve absolute paths
CONFIG="$(realpath "$CONFIG")"
PROJECT_DIR="$(realpath "$PROJECT_DIR")"

# --- Load conda for auto-detection (login node) ---
module load conda/miniforge3/24.11.3-0 2>/dev/null || true
conda activate tRNA-seq 2>/dev/null || true

# --- Auto-detect N_SAMPLES ---
if [ -n "${3:-}" ]; then
    N_SAMPLES="$3"
else
    SAMPLE_LIST=$(python3 -c "
import yaml, os
with open('$CONFIG') as f:
    cfg = yaml.safe_load(f)
sl = cfg.get('sample_list', '')
if not os.path.isabs(sl):
    sl = os.path.join('$PROJECT_DIR', sl)
print(sl)
")
    N_SAMPLES=$(python3 -c "
import pandas as pd
df = pd.read_excel('$SAMPLE_LIST')
print(len(df))
")
    echo "Auto-detected $N_SAMPLES samples from $SAMPLE_LIST"
fi

MAX_ARRAY_IDX=$((N_SAMPLES - 1))

# --- Ensure jobOutput directory exists ---
mkdir -p /home/ruv988/jobOutput

echo "=============================================================="
echo "tRNA-charge-seq SLURM Pipeline"
echo "=============================================================="
echo "Config:         $CONFIG"
echo "Project dir:    $PROJECT_DIR"
echo "Samples:        $N_SAMPLES"
echo "Max concurrent: $MAX_CONCURRENT"
echo "Job logs:       /home/ruv988/jobOutput/tseq-*"
echo "=============================================================="

# --- Job 0: Stages 0a + 0b (merge + BC split) ---
# Single job, processes all input file pairs
JOB0=$(sbatch --parsable \
    "${SCRIPT_DIR}/stage0ab.job" "$CONFIG" "$PROJECT_DIR")
echo "  JOB0 (stages 0a+0b):  $JOB0"

# --- Job 1: Stages 0c + 1 per sample (array with throttle) ---
# Each task processes one sample: UMI trim + SWIPE alignment
JOB1=$(sbatch --parsable \
    --dependency=afterok:${JOB0} \
    --array=0-${MAX_ARRAY_IDX}%${MAX_CONCURRENT} \
    "${SCRIPT_DIR}/stage0c_1.job" "$CONFIG" "$PROJECT_DIR")
echo "  JOB1 (stages 0c+1):   $JOB1  [array 0-${MAX_ARRAY_IDX}%${MAX_CONCURRENT}]"

# --- Job 2: Stages 2 + 3 + 5 (aggregation) ---
# Single job: stats collection + charge quantification + QC report
JOB2=$(sbatch --parsable \
    --dependency=afterok:${JOB1} \
    "${SCRIPT_DIR}/stage2_5.job" "$CONFIG" "$PROJECT_DIR")
echo "  JOB2 (stages 2+3+5):  $JOB2"

echo ""
echo "Pipeline submitted! Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${JOB0},${JOB1},${JOB2} --format=JobID,JobName,State,Elapsed,MaxRSS"
echo ""
echo "You will receive email notifications at russel_vincent@hms.harvard.edu"
echo "for BEGIN, END, and FAIL of each job."
