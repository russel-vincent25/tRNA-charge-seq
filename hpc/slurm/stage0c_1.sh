#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:45:00
# ==============================================================================
# Stage 0c (UMI Trimming) + Stage 1 (SWIPE Alignment) — per sample
#
# Uses SLURM_ARRAY_TASK_ID as --sample-index.
# Each array task processes ONE sample independently.
#
# Resources: 4 CPU, 8G RAM, 45 min per sample
# Bottleneck stage: SWIPE alignment (~2-3 min/sample with 4 CPUs)
# ==============================================================================

set -euo pipefail

CONFIG="$1"
PROJECT_DIR="$2"
SAMPLE_IDX="${SLURM_ARRAY_TASK_ID}"

echo "=== tRNA-charge-seq: Stage 0c + 1 (sample ${SAMPLE_IDX}) ==="
echo "Config:  $CONFIG"
echo "Project: $PROJECT_DIR"
echo "Node:    $(hostname)"
echo "CPUs:    ${SLURM_CPUS_PER_TASK:-4}"
echo "Start:   $(date)"
echo ""

# Activate conda environment
module load conda 2>/dev/null || true
source activate tRNA-seq 2>/dev/null || conda activate tRNA-seq

python -m trnaseq pipeline \
    --config "$CONFIG" \
    --project-dir "$PROJECT_DIR" \
    --n-jobs "${SLURM_CPUS_PER_TASK:-4}" \
    --stages 0c,1 \
    --sample-index "$SAMPLE_IDX"

echo ""
echo "Done: $(date)"
