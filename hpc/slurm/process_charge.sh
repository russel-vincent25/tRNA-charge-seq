#!/bin/bash

################################################################################
# SLURM Array Job Template: Process tRNA-seq Charge Data on HMS O2 Cluster
#
# This script processes charge quantification for multiple tRNA-seq samples
# in parallel using SLURM array jobs on the HMS O2 HPC cluster.
#
# Usage:
#   sbatch process_charge.sh
#   sbatch -J custom_name -o logs/charge_%a.log process_charge.sh
#
# Configuration:
#   Modify SLURM parameters below to fit your needs.
#   Array job will run for each sample ID (0-255 for PJ39).
#
################################################################################

# SLURM Configuration
#SBATCH --job-name=pj39_charge
#SBATCH --partition=short                # Use 'short' for quick jobs, 'medium' for longer
#SBATCH --time=01:00:00                  # Time limit (HH:MM:SS)
#SBATCH --mem=8G                         # Memory per task
#SBATCH --cpus-per-task=4                # CPUs per task
#SBATCH --array=0-255                    # Array job: 256 samples (0-255)
#SBATCH --output=logs/charge_%a.log      # Log file for each task
#SBATCH --error=logs/charge_%a.err       # Error file for each task
#SBATCH --mail-type=END,FAIL             # Email on completion or failure
#SBATCH --mail-user=your_email@example.com

################################################################################
# Setup and Configuration
################################################################################

# Stop on first error
set -e

# Load necessary modules (adjust based on your environment)
module load conda/23.01.0
source activate trnaseq  # Replace with your conda environment name

# Project configuration
PROJECT_DIR="/path/to/tRNA-charge-seq"              # Update this path
DATA_DIR="${PROJECT_DIR}/data"
STATS_DIR="${DATA_DIR}/stats"
CHARGE_OUTPUT_DIR="${DATA_DIR}/charge"
LOGS_DIR="${PROJECT_DIR}/logs"

# Create necessary directories
mkdir -p "${CHARGE_OUTPUT_DIR}"
mkdir -p "${LOGS_DIR}"

# Sample configuration
SAMPLE_FILE="${PROJECT_DIR}/sample_list.txt"        # File with sample names
SAMPLE_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_FILE}")

################################################################################
# Validation
################################################################################

if [ -z "${SAMPLE_ID}" ]; then
    echo "ERROR: Failed to get sample ID for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [ ! -d "${STATS_DIR}" ]; then
    echo "ERROR: Stats directory not found: ${STATS_DIR}"
    exit 1
fi

echo "=========================================="
echo "Processing Sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}"
echo "=========================================="

################################################################################
# Main Processing
################################################################################

# Use Python to process the charge data
python3 << 'PYTHON_SCRIPT'
import sys
import os
from pathlib import Path

# Add project to path
sys.path.insert(0, os.environ['PROJECT_DIR'])

import pandas as pd
from trnaseq.io import tRNAseqDataStore

# Configuration
sample_id = os.environ['SAMPLE_ID']
project_dir = Path(os.environ['PROJECT_DIR'])
data_dir = Path(os.environ['DATA_DIR'])

# Initialize storage
store = tRNAseqDataStore(data_dir)

print(f"Loading stats for sample: {sample_id}")

try:
    # Load stats for this sample
    stats_df = pd.read_csv(
        data_dir / 'stats' / f'stats_{sample_id}.csv',
        low_memory=False
    )

    print(f"Loaded {len(stats_df)} records for {sample_id}")

    # Example: Calculate charge for each tRNA
    # (Replace this with your actual charge quantification logic)
    charge_df = stats_df.groupby('tRNA_annotation').agg({
        'count': 'sum',
        'UMIcount': 'sum',
    }).reset_index()

    charge_df.rename(columns={'count': 'total_count'}, inplace=True)
    charge_df['sample_id'] = sample_id

    # Save charge data
    output_path = store.save_charge_data(charge_df, sample_id=sample_id)
    print(f"Saved charge data to: {output_path}")
    print(f"Output shape: {charge_df.shape}")

    # Verify the file was created
    if output_path.exists():
        file_size_mb = output_path.stat().st_size / (1024 * 1024)
        print(f"File size: {file_size_mb:.2f} MB")
        print("SUCCESS")
        sys.exit(0)
    else:
        print("ERROR: Output file was not created")
        sys.exit(1)

except Exception as e:
    print(f"ERROR processing {sample_id}: {str(e)}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)

PYTHON_SCRIPT

################################################################################
# Post-Processing
################################################################################

# Capture exit code
EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Task ${SLURM_ARRAY_TASK_ID} (${SAMPLE_ID}) completed successfully"
else
    echo "Task ${SLURM_ARRAY_TASK_ID} (${SAMPLE_ID}) failed with exit code ${EXIT_CODE}"
fi

exit ${EXIT_CODE}
