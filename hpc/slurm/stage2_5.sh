#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
# ==============================================================================
# Stage 2 (Stats Collection) + Stage 3 (Charge Quantification) + Stage 5 (QC)
#
# Aggregation step — processes ALL samples together.
# Reads per-sample alignment JSONs, computes stats, charge, and QC report.
#
# Resources: 16 CPU, 64G RAM, 2 hours
# For 256 samples: ~30 min stats + ~3 min charge + ~1 min QC
# ==============================================================================

set -euo pipefail

CONFIG="$1"
PROJECT_DIR="$2"

echo "=== tRNA-charge-seq: Stage 2 + 3 + 5 ==="
echo "Config:  $CONFIG"
echo "Project: $PROJECT_DIR"
echo "Node:    $(hostname)"
echo "CPUs:    ${SLURM_CPUS_PER_TASK:-16}"
echo "Start:   $(date)"
echo ""

# Activate conda environment
module load conda 2>/dev/null || true
source activate tRNA-seq 2>/dev/null || conda activate tRNA-seq

python -m trnaseq pipeline \
    --config "$CONFIG" \
    --project-dir "$PROJECT_DIR" \
    --n-jobs "${SLURM_CPUS_PER_TASK:-16}" \
    --stages 2,3,5

echo ""
echo "Done: $(date)"
