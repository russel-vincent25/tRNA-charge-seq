# tRNA-charge-seq Environment Setup Guide

This guide will help you set up the correct environment for running the tRNA-charge-seq pipeline.

## Quick Start

```bash
# 1. Create the conda environment
conda env create -f environment.yml

# 2. Activate the environment
conda activate tRNA-seq

# 3. Validate the installation
python validate_environment.py
```

If all checks pass ✅, you're ready to run the pipeline!

---

## Why This Environment File?

This `environment.yml` was created based on comprehensive beta testing (2026-02-11) that identified critical package compatibility issues:

### 🐛 Issues Resolved

**Issue #1: scipy/seaborn Incompatibility**
- **Problem:** scipy 1.15.3 breaks seaborn imports
- **Solution:** Pin scipy <1.15.0
- **Constraint:** `scipy>=1.11.0,<1.15.0`

**Issue #2: pandas/numpy Excel Reading Bug**
- **Problem:** pandas 2.3.3 + numpy 1.26.4 cannot read .xlsx files
- **Solution:** Pin pandas <2.2.0 and numpy <1.26.0
- **Constraints:**
  - `numpy>=1.24.0,<1.26.0`
  - `pandas>=2.0.0,<2.2.0`

**Issue #3: Missing Bioinformatics Tools**
- **Solution:** Include AdapterRemoval, SWIPE, BLAST, ImageMagick in environment
- All tools now installed via conda

---

## Detailed Setup Instructions

### 1. Prerequisites

- **Conda or Miniconda** installed ([Installation guide](https://docs.conda.io/en/latest/miniconda.html))
- **Disk space:** ~2-3 GB for the environment
- **Operating System:** Linux or macOS (tested on Linux)

### 2. Create Environment

```bash
# Navigate to the repository
cd /path/to/tRNA-charge-seq

# Create environment from file
conda env create -f environment.yml
```

This will:
- Install Python 3.10
- Install all required Python packages with correct versions
- Install bioinformatics tools (AdapterRemoval, SWIPE, BLAST)
- Install visualization tools (ImageMagick)

**Expected time:** 5-10 minutes

### 3. Activate Environment

```bash
conda activate tRNA-seq
```

You should see `(tRNA-seq)` prefix in your terminal prompt.

### 4. Validate Installation

Run the validation script:

```bash
python validate_environment.py
```

**Expected output:**
```
======================================================================
tRNA-charge-seq Environment Validation
======================================================================

=== Python Packages ===
✅ numpy: 1.24.4
✅ pandas: 2.0.3
✅ scipy: 1.11.4
✅ seaborn: 0.13.0
✅ matplotlib: 3.10.8
... (more packages)

=== Bioinformatics Tools ===
✅ AdapterRemoval: AdapterRemoval ver. 2.3.4
✅ SWIPE: SWIPE 2.1.1
✅ BLAST: makeblastdb: 2.5.0+
✅ ImageMagick: Version: ImageMagick 7.x

=== Testing Critical Imports ===
✅ seaborn imports successfully (scipy compatibility OK)
✅ pandas imports and works

======================================================================
✅ ALL CHECKS PASSED - Environment is ready!
======================================================================
```

### 5. Test with Example Notebook

```bash
# Start JupyterLab
jupyter lab

# Navigate to and open:
# projects/example/process_data_update.ipynb
```

Run the first few cells to confirm everything works.

---

## Troubleshooting

### ❌ Validation Script Reports Failures

**Problem:** Some packages fail version checks

**Solution:**
```bash
# Remove the environment
conda env remove -n tRNA-seq

# Recreate from scratch
conda env create -f environment.yml

# Validate again
conda activate tRNA-seq
python validate_environment.py
```

### ❌ seaborn Import Fails

**Error:** `TypeError: All ufuncs must have type 'numpy.ufunc'`

**Solution:**
```bash
conda activate tRNA-seq
pip install scipy==1.11.4 seaborn==0.13.0 --force-reinstall
```

### ❌ pandas Cannot Read Excel Files

**Error:** `TypeError: Cannot convert numpy.ndarray to numpy.ndarray`

**Solution:**
```bash
conda activate tRNA-seq
pip install numpy==1.24.4 pandas==2.0.3 --force-reinstall
```

### ❌ AdapterRemoval Not Found

**Problem:** Bioinformatics tools not in PATH

**Solution:**
```bash
conda activate tRNA-seq
conda install -c bioconda adapterremoval swipe blast
```

---

## Updating the Environment

If the `environment.yml` file is updated:

```bash
# Activate environment
conda activate tRNA-seq

# Update with new specifications
conda env update -f environment.yml --prune

# Validate
python validate_environment.py
```

---

## Alternative: Creating Environment Manually

If you prefer to create the environment manually:

```bash
# Create base environment
conda create -n tRNA-seq python=3.10 -y

# Activate
conda activate tRNA-seq

# Install packages with constraints
conda install -c conda-forge -c bioconda \
  'numpy>=1.24.0,<1.26.0' \
  'pandas>=2.0.0,<2.2.0' \
  'scipy>=1.11.0,<1.15.0' \
  seaborn=0.13.0 \
  matplotlib jupyterlab biopython \
  mpire logomaker openpyxl wand natsort \
  adapterremoval swipe blast imagemagick -y

# Validate
python validate_environment.py
```

---

## Environment Details

**Environment name:** tRNA-seq
**Python version:** 3.10
**Total packages:** ~200+ (including dependencies)
**Tested on:** Ubuntu 22.04 (WSL2), Linux 6.6.87
**Beta test date:** 2026-02-11
**Test status:** ✅ PASSED (97% success rate)

**Key metrics from beta testing:**
- AdapterRemoval merge rate: 97.32%
- BC mapping rate: 99.18% CCA+CC purity
- UMI trimming: 98.24% valid UMIs
- Alignment rate: >95%

---

## Getting Help

**If you encounter issues:**

1. **Check validation output:**
   ```bash
   python validate_environment.py
   ```

2. **Review beta test reports:**
   ```bash
   cat .claude/agent_sessions/notebook_beta_test_2026-02-11/issues_found.md
   ```

3. **Check package versions:**
   ```bash
   conda list
   ```

4. **Report issues:**
   - Include validation script output
   - Include `conda list` output
   - Describe the error message

---

## References

- **Beta Test Report:** `.claude/agent_sessions/notebook_beta_test_2026-02-11/`
- **Issue Tracker:** `.claude/agent_sessions/notebook_beta_test_2026-02-11/issues_found.md`
- **Test Summary:** `.claude/agent_sessions/notebook_beta_test_2026-02-11/test_summary.md`

---

**Last Updated:** 2026-02-11
**Status:** Production-ready ✅
