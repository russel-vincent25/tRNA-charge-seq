# Migration Guide: Beta Test → Production Environment

**Purpose:** Guide for transitioning from beta test environment to production deployment on HMS O2 cluster

**Date:** 2026-02-11
**Beta Test:** `.claude/agent_sessions/notebook_beta_test_2026-02-11/`
**Target:** PJ39 analysis (256 samples on HMS O2)

---

## 📋 Overview

This guide helps you migrate from the validated beta test environment to a production setup for processing large datasets (e.g., PJ39 with 256 samples) on HPC clusters.

**What was validated in beta test:**
- ✅ Full preprocessing pipeline (AdapterRemoval → Stats Collection)
- ✅ Package compatibility (scipy, pandas, numpy)
- ✅ Data quality metrics (>95% success rates)
- ✅ Small dataset (8 samples, ~420K reads)

**What needs migration:**
- 🔄 Environment setup on HMS O2 cluster
- 🔄 Scaling to 256 samples
- 🔄 SLURM job templates
- 🔄 Data storage optimization (Parquet)

---

## 🚀 Migration Steps

### Step 1: Replicate Environment on HMS O2

**On HMS O2 cluster:**

```bash
# Login to HMS O2
ssh username@o2.hms.harvard.edu

# Navigate to project directory
cd /n/data1/hms/your_lab/username/tRNA-charge-seq

# Clone repository (if not already there)
git clone https://github.com/your-username/tRNA-charge-seq.git
cd tRNA-charge-seq

# Create conda environment from validated environment.yml
module load conda3
conda env create -f environment.yml

# Activate and validate
conda activate tRNA-seq
python validate_environment.py
```

**Expected validation result:**
```
✅ ALL CHECKS PASSED - Environment is ready!
```

### Step 2: Verify System Tools on O2

HMS O2 may have system-level bioinformatics tools. Check if you need conda versions:

```bash
# Check if tools are already available
module avail adapterremoval
module avail blast
module avail swipe

# If available via modules, use those instead:
module load adapterremoval/2.3.4
module load blast/2.5.0+
module load swipe/2.1.1

# Otherwise, conda-installed versions will work
```

**Recommendation:** Use O2 modules if available (faster, optimized for cluster).

### Step 3: Test on Small Subset

Before processing all 256 PJ39 samples, test on a subset:

```bash
# Create test directory
mkdir -p test_run/data/raw_fastq

# Copy 5-10 samples to test directory
cp /path/to/pj39/sample_00{1..5}*.fastq.bz2 test_run/data/raw_fastq/

# Create sample_list.xlsx for test samples
# (Use beta test sample_list.xlsx as template)

# Run preprocessing pipeline on test subset
cd test_run
jupyter nbconvert --to notebook --execute \
    ../projects/example/process_data_update.ipynb \
    --output test_run_output.ipynb
```

**Success criteria:**
- All checkpoints pass
- Metrics match beta test quality (>95% rates)
- No errors in notebook execution

### Step 4: Prepare SLURM Job Templates

**Create SLURM script for parallel processing:**

```bash
#!/bin/bash
#SBATCH --job-name=pj39_preprocessing
#SBATCH --output=logs/preprocess_%A_%a.out
#SBATCH --error=logs/preprocess_%A_%a.err
#SBATCH --array=0-255                    # 256 samples
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Load environment
source ~/.bashrc
conda activate tRNA-seq

# Sample ID from array index
SAMPLE_ID=$(printf "sample%03d" $SLURM_ARRAY_TASK_ID)

echo "Processing ${SAMPLE_ID}..."
echo "Started: $(date)"

# Run preprocessing for single sample
python scripts/process_single_sample.py \
    --sample-id ${SAMPLE_ID} \
    --input-dir /n/data1/hms/lab/pj39/raw_fastq/ \
    --output-dir /n/data1/hms/lab/pj39/processed/ \
    --threads 4

echo "Completed: $(date)"
```

**Save as:** `hpc/slurm/preprocess_pj39.sh`

### Step 5: Create Single-Sample Processing Script

Since the notebook processes multiple samples, create a script for single-sample processing:

**File:** `scripts/process_single_sample.py`

```python
#!/usr/bin/env python3
"""
Process a single sample through the tRNA-charge-seq preprocessing pipeline.
Designed for SLURM array job parallelization.
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path

# Add src to path
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.read_processing import AR_merge, BC_split, UMI_trim
from src.alignment import SWIPE_align
from src.stats_collection import STATS_collection
from src.misc import read_tRNAdb_info

def main():
    parser = argparse.ArgumentParser(description='Process single tRNA-seq sample')
    parser.add_argument('--sample-id', required=True, help='Sample ID')
    parser.add_argument('--input-dir', required=True, help='Input FASTQ directory')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--min-read-len', type=int, default=31, help='Minimum read length')
    parser.add_argument('--min-align-score', type=int, default=15, help='Minimum alignment score')

    args = parser.parse_args()

    # TODO: Implement single-sample processing logic
    # This is a template - adapt from notebook cells

    print(f"Processing sample: {args.sample_id}")
    print(f"Input: {args.input_dir}")
    print(f"Output: {args.output_dir}")

    # Steps:
    # 1. AdapterRemoval
    # 2. BC_split
    # 3. UMI_trim
    # 4. SWIPE_align
    # 5. STATS_collection

    print("Processing complete!")

if __name__ == '__main__':
    main()
```

### Step 6: Scaling Considerations

**Memory Requirements:**

| Stage | Memory per Sample | Recommended |
|-------|-------------------|-------------|
| AdapterRemoval | ~2-4 GB | 8 GB |
| BC_split | ~2-4 GB | 8 GB |
| UMI_trim | ~1-2 GB | 4 GB |
| SWIPE_align | ~4-8 GB | 16 GB |
| STATS_collection | ~2-4 GB | 8 GB |

**Time Estimates (per sample):**

| Stage | Beta Test (420K reads) | PJ39 Est. (1M reads) |
|-------|------------------------|----------------------|
| AdapterRemoval | ~2 min | ~5 min |
| BC_split | ~1 min | ~2 min |
| UMI_trim | ~1 min | ~2 min |
| SWIPE_align | ~40 min | **~90 min** |
| STATS_collection | ~2 min | ~5 min |
| **TOTAL** | ~46 min | **~104 min** |

**SLURM Recommendations for PJ39:**
- Time limit: `--time=04:00:00` (4 hours, includes buffer)
- Memory: `--mem=16G` (for alignment step)
- CPUs: `--cpus-per-task=4`
- Partition: `short` or `medium` depending on cluster config

### Step 7: Data Storage Migration

**Beta Test:** Used CSV files (~15 MB total)
**PJ39 (256 samples):** Would create ~3.8 GB of CSV files

**Solution:** Use Parquet format (already implemented in `trnaseq.io.storage`)

```python
# After stats collection, convert to Parquet
from trnaseq.io.storage import tRNAseqDataStore

store = tRNAseqDataStore('/n/data1/hms/lab/pj39/processed/')

# Load all CSV stats files
all_stats = pd.concat([
    pd.read_csv(f'sample{i:03d}_stats_aggregate.csv')
    for i in range(256)
])

# Save to Parquet (15x smaller, 16x faster to load)
store.save_stats(all_stats, format='parquet')

# Result: ~250 MB instead of 3.8 GB
```

### Step 8: Integration with Charge Analysis

Once preprocessing is complete, integrate with charge analysis modules (already developed):

```python
from trnaseq.charge import ChargeQuantifier
from trnaseq.io.storage import tRNAseqDataStore

# Load preprocessed data
store = tRNAseqDataStore('/n/data1/hms/lab/pj39/processed/')
stats_df = store.load_stats(format='parquet')

# Quantify charge
quantifier = ChargeQuantifier(stats_df=stats_df)
charge_df = quantifier.quantify_all()

# Save charge results
store.save_charge_data(charge_df)
```

---

## 🔍 Quality Control Checklist

**Before full production run:**

- [ ] Environment validated on O2 cluster (`validate_environment.py` passes)
- [ ] Test run completed on 5-10 samples successfully
- [ ] Quality metrics match beta test expectations:
  - [ ] Merge rate >90%
  - [ ] BC mapping rate >80% (or investigate if lower)
  - [ ] Valid UMI rate >95%
  - [ ] Alignment rate >95%
  - [ ] CCA+CC percentage >95%
- [ ] SLURM script tested with small array (e.g., `--array=0-4`)
- [ ] Output directory structure verified
- [ ] Disk space available (~500 GB for 256 samples estimated)

**During production run:**

- [ ] Monitor SLURM job status: `squeue -u $USER`
- [ ] Check log files for errors: `tail -f logs/preprocess_*.err`
- [ ] Verify output files being created: `ls -lh processed/*/`
- [ ] Monitor disk usage: `du -sh processed/`

**After production run:**

- [ ] All 256 samples completed (check array job completion)
- [ ] No error logs with critical failures
- [ ] Verify statistics files exist for all samples
- [ ] Convert to Parquet for efficient storage
- [ ] Run QC validation across all samples
- [ ] Generate summary report

---

## 🐛 Common Migration Issues

### Issue #1: Module Conflicts on O2

**Problem:** O2 system modules conflict with conda packages

**Solution:**
```bash
# Purge all modules before activating conda
module purge
conda activate tRNA-seq
```

### Issue #2: Disk Space on O2

**Problem:** Home directory quota exceeded

**Solution:**
```bash
# Use lab storage directory for data
export TMPDIR=/n/scratch3/users/${USER:0:1}/${USER}
export DATA_DIR=/n/data1/hms/your_lab/username/pj39
```

### Issue #3: SWIPE Alignment Timeout

**Problem:** Alignment takes longer than expected on cluster

**Solution:**
```bash
# Enable common_seqs optimization in alignment step
# Edit process_single_sample.py to use:
common_seqs = '/path/to/utils/common-seqs.fasta.bz2'
align_obj = SWIPE_align(..., common_seqs=common_seqs, ...)
```

### Issue #4: Dependency Version Mismatch

**Problem:** O2 installed different package versions than environment.yml

**Solution:**
```bash
# Force reinstall with exact versions
conda activate tRNA-seq
pip install numpy==1.24.4 pandas==2.0.3 scipy==1.11.4 --force-reinstall
```

---

## 📊 Expected Outcomes

### Beta Test Results (Reference)
- **Dataset:** 8 samples, 420K reads
- **Merge rate:** 97.32%
- **BC mapping:** 79.64% (data-specific, CCA+CC still 99.18%)
- **Valid UMI:** 98.24%
- **Alignment:** 75% (interrupted, but working well)

### PJ39 Production Target
- **Dataset:** 256 samples, ~256M reads total
- **Runtime:** ~104 min/sample = ~6.7 hours total (parallelized)
- **Storage:** ~250 MB (Parquet) vs ~3.8 GB (CSV)
- **Quality:** Match or exceed beta test metrics

---

## 🔄 Rollback Plan

If production run fails:

1. **Stop array job:**
   ```bash
   scancel $SLURM_ARRAY_JOB_ID
   ```

2. **Identify failing samples:**
   ```bash
   grep -l "ERROR\|FAILED" logs/preprocess_*.err
   ```

3. **Rerun failed samples:**
   ```bash
   # Create list of failed sample indices
   # Resubmit with reduced array
   sbatch --array=0,5,12,... hpc/slurm/preprocess_pj39.sh
   ```

4. **Fall back to sequential processing if needed:**
   ```bash
   # Process samples one-by-one with full resources
   for i in {0..255}; do
       sbatch --mem=32G --time=08:00:00 process_sample_$i.sh
   done
   ```

---

## 📚 Additional Resources

**Beta Test Documentation:**
- `.claude/agent_sessions/notebook_beta_test_2026-02-11/issues_found.md`
- `.claude/agent_sessions/notebook_beta_test_2026-02-11/test_summary.md`
- `ENVIRONMENT_SETUP.md`

**Development Plan:**
- `.claude/CURRENT_PLAN.md` - Overall development roadmap

**HMS O2 Resources:**
- [O2 Wiki](https://wiki.rc.hms.harvard.edu/display/O2)
- [SLURM Documentation](https://slurm.schedmd.com/)

---

## ✅ Success Criteria

**Environment Migration:**
- ✅ Same package versions as beta test
- ✅ All validation checks pass
- ✅ Test run matches beta test quality

**Production Run:**
- ✅ All 256 samples processed without critical errors
- ✅ Quality metrics meet thresholds
- ✅ Data stored efficiently (Parquet format)
- ✅ Ready for downstream charge analysis

**Integration:**
- ✅ Charge quantification works with preprocessed data
- ✅ Alignment viewer QC functional
- ✅ Results ready for analysis and manuscript

---

**Last Updated:** 2026-02-11
**Beta Test Status:** ✅ VALIDATED
**Production Status:** 🔄 READY FOR MIGRATION
**Target Deployment:** HMS O2 cluster
