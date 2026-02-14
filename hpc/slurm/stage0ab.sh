#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# ==============================================================================
# Stage 0a (Adapter Removal + Merging) + Stage 0b (Barcode Splitting)
#
# Processes ALL input file pairs — not per-sample.
# Typically fast: ~5 min for 8 samples, ~35 min for 256 samples.
#
# Resources: 4 CPU, 16G RAM, 1 hour
# ==============================================================================

set -euo pipefail

CONFIG="$1"
PROJECT_DIR="$2"

echo "=== tRNA-charge-seq: Stage 0a + 0b ==="
echo "Config:  $CONFIG"
echo "Project: $PROJECT_DIR"
echo "Node:    $(hostname)"
echo "CPUs:    ${SLURM_CPUS_PER_TASK:-4}"
echo "Start:   $(date)"
echo ""

# Activate conda environment
# On HMS O2: module load conda, then activate
module load conda 2>/dev/null || true
source activate tRNA-seq 2>/dev/null || conda activate tRNA-seq

python -m trnaseq pipeline \
    --config "$CONFIG" \
    --project-dir "$PROJECT_DIR" \
    --n-jobs "${SLURM_CPUS_PER_TASK:-4}" \
    --stages 0a,0b

echo ""
echo "Done: $(date)"
