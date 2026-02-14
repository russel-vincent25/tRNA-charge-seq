#!/bin/bash
# ==============================================================================
# Master launcher for tRNA-charge-seq pipeline on HMS O2 (SLURM)
# ==============================================================================
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
# Job chain:
#   JOB0 (stage0ab)  →  JOB1 (stage0c_1, array)  →  JOB2 (stage2_5)
#
# Example:
#   bash submit_pipeline.sh /path/to/config.yaml /path/to/PJ39/ 256 32
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

# --- Create logs directory ---
LOGS_DIR="${PROJECT_DIR}/logs"
mkdir -p "$LOGS_DIR"

echo "=============================================================="
echo "tRNA-charge-seq SLURM Pipeline"
echo "=============================================================="
echo "Config:         $CONFIG"
echo "Project dir:    $PROJECT_DIR"
echo "Samples:        $N_SAMPLES"
echo "Max concurrent: $MAX_CONCURRENT"
echo "Logs:           $LOGS_DIR"
echo "=============================================================="

# --- Job 0: Stages 0a + 0b (merge + BC split) ---
# Single job, processes all input file pairs
JOB0=$(sbatch --parsable \
    --job-name=tseq-0ab \
    --output="${LOGS_DIR}/stage0ab_%j.out" \
    --error="${LOGS_DIR}/stage0ab_%j.err" \
    "${SCRIPT_DIR}/stage0ab.sh" "$CONFIG" "$PROJECT_DIR")
echo "  JOB0 (stages 0a+0b):  $JOB0"

# --- Job 1: Stages 0c + 1 per sample (array with throttle) ---
# Each task processes one sample: UMI trim + SWIPE alignment
JOB1=$(sbatch --parsable \
    --dependency=afterok:${JOB0} \
    --job-name=tseq-0c1 \
    --array=0-${MAX_ARRAY_IDX}%${MAX_CONCURRENT} \
    --output="${LOGS_DIR}/stage0c1_%A_%a.out" \
    --error="${LOGS_DIR}/stage0c1_%A_%a.err" \
    "${SCRIPT_DIR}/stage0c_1.sh" "$CONFIG" "$PROJECT_DIR")
echo "  JOB1 (stages 0c+1):   $JOB1  [array 0-${MAX_ARRAY_IDX}%${MAX_CONCURRENT}]"

# --- Job 2: Stages 2 + 3 + 5 (aggregation) ---
# Single job: stats collection + charge quantification + QC report
JOB2=$(sbatch --parsable \
    --dependency=afterok:${JOB1} \
    --job-name=tseq-235 \
    --output="${LOGS_DIR}/stage235_%j.out" \
    --error="${LOGS_DIR}/stage235_%j.err" \
    "${SCRIPT_DIR}/stage2_5.sh" "$CONFIG" "$PROJECT_DIR")
echo "  JOB2 (stages 2+3+5):  $JOB2"

echo ""
echo "Pipeline submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${JOB0},${JOB1},${JOB2} --format=JobID,State,Elapsed,MaxRSS"
echo "  tail -f ${LOGS_DIR}/*.out"
