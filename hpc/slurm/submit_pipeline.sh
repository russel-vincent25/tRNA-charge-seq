#!/bin/bash
# ==============================================================================
# Master launcher for tRNA-charge-seq pipeline on HMS O2 (SLURM)
# ==============================================================================
#
# This script runs on the LOGIN NODE. It submits 3 sbatch jobs that execute
# on compute nodes with dependency chaining:
#
#   JOB0 (stage0ab.job)  →  JOB1 (stage0c_1.job, array)  →  JOB2 (stage2_6.job)
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

# --- Resource scaling based on N_SAMPLES ---
THREADS_PER_JOB=2

# Stage 0ab: process all file pairs; scale CPUs with sample count, cap at 16
# RAM: AR_merge + BC_split are lightweight stream processors (~500MB for 24 samples)
CPUS_0AB=$(( N_SAMPLES < 16 ? (N_SAMPLES > 2 ? N_SAMPLES : 2) : 16 ))
N_JOBS_0AB=$(( CPUS_0AB / THREADS_PER_JOB ))
# ~6 min per batch (AR merge + BC split), minimum 2h, cap 12h
TIME_0AB_MIN=$(( (N_SAMPLES / N_JOBS_0AB + 1) * 6 ))
if [ $TIME_0AB_MIN -lt 120 ]; then TIME_0AB_MIN=120; fi
if [ $TIME_0AB_MIN -gt 720 ]; then TIME_0AB_MIN=720; fi
TIME_0AB=$(printf "%02d:%02d:00" $((TIME_0AB_MIN / 60)) $((TIME_0AB_MIN % 60)))
MEM_0AB="4G"

# Stage 0c+1: per-sample array; SWIPE is I/O-bound (~50% CPU efficiency)
# 2 CPUs sufficient; RAM peaks at ~1GB even for large samples
CPUS_0C1=2
TIME_0C1="08:00:00"
MEM_0C1="2G"

# Stage 2+3+5+6+7: aggregation; I/O-bound (reading JSONs, writing CSVs)
# 42 min for 24 samples with 4 CPUs → scales ~linearly with N_SAMPLES/CPUS
# 16 CPUs keeps 264 samples under 4h; RAM scales with sample count
CPUS_2356=16
if [ $N_SAMPLES -lt 64 ]; then MEM_2356="16G"; else MEM_2356="32G"; fi
# ~1.75 min per sample per CPU → estimate wall time, minimum 1h, cap 12h
TIME_2356_MIN=$(( (N_SAMPLES * 2 / CPUS_2356 + 1) * 3 ))
if [ $TIME_2356_MIN -lt 60 ]; then TIME_2356_MIN=60; fi
if [ $TIME_2356_MIN -gt 720 ]; then TIME_2356_MIN=720; fi
TIME_2356=$(printf "%02d:%02d:00" $((TIME_2356_MIN / 60)) $((TIME_2356_MIN % 60)))

echo "=============================================================="
echo "tRNA-charge-seq SLURM Pipeline"
echo "=============================================================="
echo "Config:         $CONFIG"
echo "Project dir:    $PROJECT_DIR"
echo "Samples:        $N_SAMPLES"
echo "Max concurrent: $MAX_CONCURRENT"
echo "Threads/job:    $THREADS_PER_JOB"
echo "Job logs:       /home/ruv988/jobOutput/tseq-*"
echo ""
echo "Resource plan:"
echo "  Stage 0ab:  ${CPUS_0AB} CPUs, ${MEM_0AB} mem, ${TIME_0AB} time (${N_JOBS_0AB} jobs × ${THREADS_PER_JOB} threads)"
echo "  Stage 0c+1: ${CPUS_0C1} CPUs, ${MEM_0C1} mem, ${TIME_0C1} time (1 job × ${CPUS_0C1} threads) × ${N_SAMPLES} tasks"
echo "  Stage 2+3+5+6+7: ${CPUS_2356} CPUs, ${MEM_2356} mem, ${TIME_2356} time"
echo "=============================================================="

# --- Job 0: Stages 0a + 0b (merge + BC split) ---
# Single job, processes all input file pairs
JOB0=$(sbatch --parsable \
    --cpus-per-task=$CPUS_0AB --mem=${MEM_0AB} -t ${TIME_0AB} \
    "${SCRIPT_DIR}/stage0ab.job" "$CONFIG" "$PROJECT_DIR" "$THREADS_PER_JOB")
echo "  JOB0 (stages 0a+0b):  $JOB0"

# --- Job 1: Stages 0c + 1 per sample (array with throttle) ---
# Each task processes one sample: UMI trim + SWIPE alignment
JOB1=$(sbatch --parsable \
    --dependency=afterok:${JOB0} \
    --array=0-${MAX_ARRAY_IDX}%${MAX_CONCURRENT} \
    --cpus-per-task=$CPUS_0C1 --mem=${MEM_0C1} -t ${TIME_0C1} \
    "${SCRIPT_DIR}/stage0c_1.job" "$CONFIG" "$PROJECT_DIR")
echo "  JOB1 (stages 0c+1):   $JOB1  [array 0-${MAX_ARRAY_IDX}%${MAX_CONCURRENT}]"

# --- Job 2: Stages 2 + 3 + 5 + 6 + 7 (aggregation) ---
# Single job: stats + charge + fragments + QC + modifications + abundance
JOB2=$(sbatch --parsable \
    --dependency=afterok:${JOB1} \
    --cpus-per-task=$CPUS_2356 --mem=${MEM_2356} -t ${TIME_2356} \
    "${SCRIPT_DIR}/stage2_6.job" "$CONFIG" "$PROJECT_DIR")
echo "  JOB2 (stages 2+3+5+6+7):  $JOB2"

echo ""
echo "Pipeline submitted! Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${JOB0},${JOB1},${JOB2} --format=JobID,JobName,State,Elapsed,MaxRSS"
echo ""
echo "You will receive email notifications at russel_vincent@hms.harvard.edu"
echo "for BEGIN, END, and FAIL of each job."
