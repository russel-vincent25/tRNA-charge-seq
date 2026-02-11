# tRNA-charge-seq Scripts

This directory contains standalone scripts for running the tRNA-charge-seq pipeline.

## Main Scripts

### `run_preprocessing.py` - Unified Preprocessing Pipeline

Combines Stages 0-2 into a single command:
- Stage 0a: Adapter Removal & Merging
- Stage 0b: Barcode Splitting
- Stage 0c: UMI Trimming
- Stage 1: SWIPE Alignment
- Stage 2: Stats Collection

**Input:** Raw FASTQ files + configuration file

**Output:**
- `inp_file_df.xlsx` - Input file summary with QC metrics
- `sample_df.xlsx` - Sample information with QC metrics
- `ALL_stats_aggregate.csv` - Combined statistics for charge quantification

**Usage:**
```bash
# Basic usage
python scripts/run_preprocessing.py \
    --config config.yaml \
    --output-dir output/ \
    --n-jobs 8

# With custom configuration
python scripts/run_preprocessing.py \
    --config projects/PJ39/config_pj39.yaml \
    --output-dir /n/data1/lab/PJ39-output/ \
    --n-jobs 16
```

**Configuration File:**
See `config_example.yaml` in repository root for a template.

**SLURM Example:**
```bash
#!/bin/bash
#SBATCH --job-name=preprocessing
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=0-04:00

module load conda3/latest
conda activate trnaseq

python scripts/run_preprocessing.py \
    --config config_pj39.yaml \
    --output-dir output/ \
    --n-jobs 16
```

**Output Structure:**
```
output/
├── preprocessing.log                # Detailed log file
├── inp_file_df.xlsx                # Input file summary
├── sample_df.xlsx                  # Sample QC metrics
├── ALL_stats_aggregate.csv         # Ready for charge quantification
└── data/
    ├── AdapterRemoval/
    ├── BC_split/
    ├── UMI_trimmed/
    ├── SWalign/
    └── stats_collection/
```

**Next Step After Preprocessing:**
```bash
# Run charge quantification (Stage 3)
python -m trnaseq.cli.commands.quantify \
    output/ALL_stats_aggregate.csv \
    --output charge_results.csv
```

---

## Configuration Files

Configuration files are in YAML format and define all pipeline parameters.

**Key sections:**
- `sample_list` - Path to sample_list.xlsx
- `index_list` - Path to index_list.xlsx with barcode definitions
- `tRNA_database` - Paths to species-specific tRNA reference databases
- `SWIPE_score_mat` - Alignment scoring matrix
- `min_read_len` - Minimum read length after merging (default: 39)
- `downsample_percentile` - Optional: even out samples
- `downsample_absolute` - Optional: cap reads per sample

**Example:**
```yaml
sample_list: "projects/PJ39/sample_list.xlsx"
index_list: "utils/index_list.xlsx"
tRNA_database:
  human: "tRNA_database/human/hg38-tRNAs.fa"
min_read_len: 39
downsample_percentile: null  # No downsampling
```

---

## Quality Control

The pipeline automatically logs QC metrics:

**Stage 0a (Adapter Removal):**
- Expected: >90% reads successfully merged

**Stage 0b (Barcode Split):**
- Expected: >90% reads mapped to a barcode
- Warning if <80%: Check for missing samples or bad adapters

**Stage 0c (UMI Trim):**
- Expected: >95% valid UMI sequences
- Expected: >90% observed vs expected UMIs
- Warning if <80%: Library prep bottleneck

**Stage 1 (Alignment):**
- Expected: >95% reads aligned
- Expected: >70% uniquely mapped
- Warning if <80%: Check reference database

**All metrics are saved in `sample_df.xlsx` for review.**

---

## Troubleshooting

### Low Barcode Mapping (<80%)

The pipeline will warn if barcode mapping is low. Check:
1. `sample_list.xlsx` includes all barcodes present in data
2. Barcode sequences in `index_list.xlsx` are correct
3. No adapter synthesis errors

### Low UMI Diversity (<80% obs/exp)

This suggests library prep bottleneck:
1. Starting material was too low (same molecules sequenced multiple times)
2. PCR over-amplification
3. Solution: Use UMI counts instead of read counts for quantification

### Low Mapping Rate (<80%)

Check:
1. Correct tRNA database for species
2. Reference database is complete
3. Sequencing quality is good

---

## Advanced: Partial Runs

To run only specific stages, modify the script or use the individual
classes from `src/`:

```python
from src.read_processing import AR_merge, BC_split, UMI_trim
from src.alignment import SWIPE_align
from src.stats_collection import STATS_collection

# Run only alignment
align_obj = SWIPE_align(dir_dict, tRNA_database, sample_df, ...)
sample_df = align_obj.run_parallel(n_jobs=8)
```

---

**For detailed examples, see:** `projects/example/process_data.ipynb`
