# tRNA-seq Analysis Package - Architecture Design
**Date:** 2026-01-23 (Updated: 2026-02-10)
**Goal:** Professional bioinformatics tool (minimap2/nextflow style) with modular design

---

## ⚠️ **Important Update (2026-02-10)**

The package structure now reflects the **complete pipeline** including preprocessing steps:
1. **Adapter Removal & Merging** (AdapterRemoval wrapper)
2. **Barcode Splitting** (Demultiplexing by adapter barcodes)
3. **UMI Trimming** (10nt UMI extraction)
4. **Alignment** (SWIPE to tRNA database)
5. **Stats Collection** (Alignment → CSV)
6. **Charge Quantification** (CSV → Charge data) ← **NEW** (implemented 2026-02-10)

See `projects/example/process_data.ipynb` for complete workflow example.

---

## Design Philosophy

### Inspiration: minimap2 + nf-core

**minimap2 principles:**
- ✅ Single executable with multiple modes
- ✅ Clear subcommands (preprocess, align, quantify, analyze)
- ✅ Can be used as CLI tool OR imported as library
- ✅ Optimized for HPC (parallel processing, memory efficient)
- ✅ Professional documentation

**nf-core principles:**
- ✅ Nextflow-ready modules
- ✅ Containerized (Docker/Singularity)
- ✅ Reproducible workflows
- ✅ Standard parameters and configs

---

## Package Structure

```
tRNA-charge-seq/                    # Main repository (your forked repo)
├── pyproject.toml                  # Modern Python packaging (PEP 517/518)
├── setup.py                        # Backwards compatibility
├── README.md
├── LICENSE
├── CHANGELOG.md
│
├── trnaseq/                        # Main package
│   ├── __init__.py
│   ├── __version__.py
│   │
│   ├── core/                       # Core functionality (current: src/)
│   │   ├── __init__.py
│   │   ├── alignment.py            # SWIPE wrapper (src/alignment.py)
│   │   ├── read_processing.py      # AR_merge, BC_split, UMI_trim (src/read_processing.py)
│   │   ├── stats_collection.py     # STATS_collection class (src/stats_collection.py)
│   │   ├── curve_fitting.py        # Decay kinetics (src/curve_fitting.py)
│   │   └── misc.py                 # Helper functions (src/misc.py)
│   │
│   ├── modifications/              # RT signature analysis
│   │   ├── __init__.py
│   │   ├── rt_signatures.py        # Extract from transcript_mutations.py
│   │   ├── modification_caller.py  # m1A, m3C, Ψ, m7G annotation
│   │   ├── modomics.py             # MODOMICS database integration
│   │   └── visualize.py            # Modification plots
│   │
│   ├── fragments/                  # Fragment analysis (from enhanced branch)
│   │   ├── __init__.py
│   │   ├── classifier.py           # 7+ fragment categories
│   │   ├── quantifier.py           # Count matrices, ratios
│   │   ├── differential.py         # Statistical tests
│   │   └── visualize.py            # Fragment plots
│   │
│   ├── charge/                     # Charge state analysis ✅ IMPLEMENTED (2026-02-10)
│   │   ├── __init__.py
│   │   ├── quantifier.py           # ChargeQuantifier class (extracts from CSV)
│   │   ├── example_usage.py        # Usage examples
│   │   └── decay_fitting.py        # Curve fitting (from src/curve_fitting.py)
│   │
│   ├── preprocessing/              # Read preprocessing (future: refactor from core/)
│   │   ├── __init__.py
│   │   ├── adapter_removal.py      # AdapterRemoval wrapper (AR_merge class)
│   │   ├── demultiplex.py          # Barcode splitting (BC_split class)
│   │   ├── umi_extraction.py       # UMI trimming (UMI_trim class)
│   │   └── quality_control.py      # QC metrics (Kmer_analysis, BC_analysis)
│   │
│   ├── visualization/              # Visualization tools ✅ IMPLEMENTED (2026-02-10)
│   │   ├── __init__.py
│   │   ├── alignment_viewer.py     # AlignmentViewer class (JSON.bz2 → plots/HTML)
│   │   ├── example_usage.py        # Usage examples
│   │   ├── README.md               # Documentation
│   │   └── (future: coverage.py, logos.py, heatmaps.py from plotting.py)
│   │
│   ├── utils/                      # Shared utilities
│   │   ├── __init__.py
│   │   ├── io.py                   # File I/O (BAM, FASTQ, CSV)
│   │   ├── parallel.py             # Parallel processing helpers
│   │   ├── validators.py           # Input validation
│   │   └── config.py               # Config file handling
│   │
│   └── cli/                        # Command-line interface
│       ├── __init__.py
│       ├── main.py                 # Entry point: trnaseq command
│       └── commands/
│           ├── __init__.py
│           ├── preprocess.py       # trnaseq preprocess
│           ├── align.py            # trnaseq align
│           ├── quantify.py         # trnaseq quantify
│           ├── modifications.py    # trnaseq modifications
│           ├── fragments.py        # trnaseq fragments
│           └── run.py              # trnaseq run (full pipeline)
│
├── workflows/                      # Nextflow/Snakemake workflows
│   ├── nextflow/
│   │   ├── main.nf
│   │   ├── nextflow.config
│   │   ├── modules/
│   │   │   ├── preprocess.nf
│   │   │   ├── align.nf
│   │   │   ├── quantify.nf
│   │   │   ├── modifications.nf
│   │   │   └── fragments.nf
│   │   └── profiles/
│   │       ├── standard.config
│   │       ├── slurm.config
│   │       └── aws.config
│   │
│   └── snakemake/
│       ├── Snakefile
│       └── rules/
│
├── containers/                     # Containerization
│   ├── Dockerfile
│   ├── Singularity.def
│   └── conda_env.yaml
│
├── configs/                        # Configuration examples
│   ├── default.yaml
│   ├── hms_o2.yaml                 # HMS O2 cluster
│   └── examples/
│       ├── modification_analysis.yaml
│       ├── fragment_analysis.yaml
│       └── full_pipeline.yaml
│
├── data/                          # Reference data
│   ├── modomics/
│   │   ├── human_modifications.json
│   │   ├── mouse_modifications.json
│   │   └── ecoli_modifications.json
│   └── examples/
│       └── test_data/
│
├── tests/                         # Test suite
│   ├── unit/
│   │   ├── test_modifications.py
│   │   ├── test_fragments.py
│   │   └── test_charge.py
│   ├── integration/
│   │   └── test_full_pipeline.py
│   └── data/
│       └── test_samples/
│
├── docs/                          # Documentation (Sphinx)
│   ├── conf.py
│   ├── index.rst
│   ├── installation.md
│   ├── quickstart.md
│   ├── user_guide/
│   │   ├── preprocessing.md
│   │   ├── alignment.md
│   │   ├── modifications.md
│   │   ├── fragments.md
│   │   └── workflows.md
│   ├── api_reference/
│   └── tutorials/
│       ├── notebook_usage.ipynb
│       └── cli_usage.md
│
├── notebooks/                     # Example notebooks (for users)
│   ├── 01_preprocessing.ipynb
│   ├── 02_alignment.ipynb
│   ├── 03_charge_analysis.ipynb
│   ├── 04_modification_analysis.ipynb
│   ├── 05_fragment_analysis.ipynb
│   └── 06_full_workflow.ipynb
│
└── scripts/                       # Helper scripts
    ├── download_modomics.py
    ├── generate_test_data.py
    └── benchmark.sh
```

---

## Dual Interface Design

### 1. Python Library (for notebooks and custom scripts)

**Import style:**
```python
# High-level API (simplified)
from trnaseq import ModificationAnalyzer, FragmentAnalyzer

# Modification analysis in notebook
mod = ModificationAnalyzer(
    bam_files=['sample1.bam', 'sample2.bam'],
    reference='hg38-tRNAs.fa',
    organism='human'
)

# Run analysis
signatures = mod.detect_rt_signatures(min_coverage=50)
modifications = mod.call_modifications(signatures, confidence=0.7)
mod.plot_heatmap(modifications, output='mod_heatmap.pdf')

# Fragment analysis
frag = FragmentAnalyzer(
    bam_files=['sample1.bam'],
    reference='hg38-tRNAs.fa'
)

fragments = frag.classify_fragments()
counts = frag.build_count_matrix(fragments)
frag.plot_size_distribution(fragments, output='frag_dist.pdf')
```

**Low-level API (for advanced users):**
```python
# Direct module access
from trnaseq.modifications import RTSignatureDetector, ModificationCaller
from trnaseq.modifications.modomics import MODOMICSDatabase

# Fine-grained control
detector = RTSignatureDetector(min_coverage=50, mismatch_threshold=0.10)
signatures = detector.analyze_bam('sample1.bam', 'hg38-tRNAs.fa')

# Custom modification calling
modomics = MODOMICSDatabase(organism='human')
caller = ModificationCaller(modomics_db=modomics)
mods = caller.call_from_signatures(signatures, confidence=0.7)
```

### 2. Command-Line Interface (for HPC and automation)

**Main command structure:**
```bash
trnaseq [OPTIONS] COMMAND [ARGS]
```

**Subcommands:**

```bash
# Preprocessing
trnaseq preprocess \
    --fastq-dir data/raw/ \
    --output-dir data/preprocessed/ \
    --adapter-r1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --threads 8

# Alignment
trnaseq align \
    --input-dir data/preprocessed/ \
    --reference refs/hg38-tRNAs.fa \
    --output-dir data/aligned/ \
    --threads 16

# Charge quantification
trnaseq quantify \
    --bam-dir data/aligned/ \
    --charge-type canonical \
    --output-dir data/quantified/

# Modification analysis
trnaseq modifications \
    --bam data/aligned/*.bam \
    --reference refs/hg38-tRNAs.fa \
    --output data/modifications/ \
    --organism human \
    --min-coverage 50 \
    --call-modifications

# Fragment analysis
trnaseq fragments classify \
    --bam data/aligned/*.bam \
    --reference refs/hg38-tRNAs.fa \
    --output data/fragments/

trnaseq fragments quantify \
    --bam data/aligned/*.bam \
    --sample-info samples.csv \
    --output data/fragments/counts/

trnaseq fragments differential \
    --counts data/fragments/counts/fragment_counts.csv \
    --sample-info samples.csv \
    --comparison Treated,Control

# Full pipeline from config
trnaseq run \
    --config configs/full_pipeline.yaml \
    --output-dir results/ \
    --log-file pipeline.log

# Background execution
trnaseq run \
    --config configs/full_pipeline.yaml \
    --background \
    --slurm \
    --account yankner_lab
```

---

## Configuration File Format

**File:** `configs/modification_analysis.yaml`

```yaml
# tRNA-seq Modification Analysis Configuration

project:
  name: RVHMS09_AD_tRNAseq
  description: Alzheimer's Disease tRNA modifications
  output_dir: results/rvhms09/

input:
  bam_dir: data/aligned/
  reference: refs/hg38-tRNAs.fa
  sample_info: metadata/samples.csv

parameters:
  organism: human
  min_coverage: 50
  mismatch_threshold: 0.10
  rt_stop_threshold: 0.20

modifications:
  call_modifications: true
  confidence_threshold: 0.7
  modomics_validation: true
  known_modifications:
    - m1A
    - m3C
    - pseudouridine
    - m7G
    - m5C

output:
  signatures_csv: rt_signatures.csv
  modifications_csv: modifications_called.csv
  heatmap_pdf: modification_heatmap.pdf
  coverage_plots: true
  logo_plots: true

resources:
  threads: 8
  memory: 16G
```

---

## Installation & Distribution

### PyPI Installation

```bash
# Stable release
pip install trnaseq-charge

# Development version
pip install git+https://github.com/russel-vincent25/tRNA-charge-seq.git

# With optional dependencies
pip install trnaseq-charge[nextflow,plotting,ml]
```

### Conda Installation

```bash
conda install -c bioconda trnaseq-charge
```

### Docker Container

```bash
# Pull from Docker Hub
docker pull russelv/trnaseq-charge:latest

# Run analysis in container
docker run -v $(pwd):/data russelv/trnaseq-charge:latest \
    trnaseq modifications \
    --bam /data/aligned/*.bam \
    --reference /data/refs/hg38-tRNAs.fa \
    --output /data/modifications/
```

### Singularity (for HPC)

```bash
# Build container
singularity build trnaseq.sif docker://russelv/trnaseq-charge:latest

# Run on cluster
singularity exec trnaseq.sif trnaseq modifications ...
```

---

## Nextflow Workflow Integration

### Modular Design

Each analysis step is a separate Nextflow module:

**File:** `workflows/nextflow/modules/modifications.nf`

```groovy
process MODIFICATIONS_DETECT {
    tag "$sample_id"
    label 'process_medium'
    conda 'bioconda::trnaseq-charge'
    container 'russelv/trnaseq-charge:latest'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_rt_signatures.csv"), emit: signatures
    path "${sample_id}_signatures.log", emit: log

    script:
    """
    trnaseq modifications \
        --bam ${bam} \
        --reference ${reference} \
        --output . \
        --min-coverage 50 \
        --sample-name ${sample_id} \
        --signatures-only
    """
}

process MODIFICATIONS_CALL {
    tag "$sample_id"
    label 'process_low'
    conda 'bioconda::trnaseq-charge'
    container 'russelv/trnaseq-charge:latest'

    input:
    tuple val(sample_id), path(signatures)
    val organism

    output:
    tuple val(sample_id), path("${sample_id}_modifications.csv"), emit: modifications
    path "${sample_id}_modifications.pdf", emit: plots

    script:
    """
    trnaseq modifications call \
        --signatures ${signatures} \
        --organism ${organism} \
        --confidence 0.7 \
        --output . \
        --sample-name ${sample_id}
    """
}
```

**Main workflow:** `workflows/nextflow/main.nf`

```groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { MODIFICATIONS_DETECT; MODIFICATIONS_CALL } from './modules/modifications'
include { FRAGMENTS_CLASSIFY; FRAGMENTS_QUANTIFY } from './modules/fragments'

workflow {
    // Input channels
    bam_ch = Channel
        .fromFilePairs(params.bam_pattern, size: 2)
        .map { id, files -> [id, files[0], files[1]] }

    reference = file(params.reference)
    organism = params.organism

    // Modification analysis
    MODIFICATIONS_DETECT(bam_ch, reference)
    MODIFICATIONS_CALL(MODIFICATIONS_DETECT.out.signatures, organism)

    // Fragment analysis
    FRAGMENTS_CLASSIFY(bam_ch, reference)
    FRAGMENTS_QUANTIFY(FRAGMENTS_CLASSIFY.out.fragments, params.sample_info)

    // Emit results
    modifications_ch = MODIFICATIONS_CALL.out.modifications
    fragments_ch = FRAGMENTS_QUANTIFY.out.counts
}
```

---

## Development Workflow

### Phase 1: Extract and Modularize (Week 1)

**Goal:** Create standalone modules from existing code

**Tasks:**
1. **Day 1-2:** Extract RT signature analysis
   - Create `trnaseq/modifications/rt_signatures.py`
   - Extract from `transcript_mutations.py::TM_analysis`
   - Keep PSCM, mutation frequency, RT stops logic
   - Add unit tests

2. **Day 3-4:** Create modification caller
   - Create `trnaseq/modifications/modification_caller.py`
   - Implement known modification profiles (m1A, m3C, Ψ, m7G)
   - Add statistical testing (binomial)
   - Download MODOMICS data

3. **Day 5:** Build high-level API
   - Create `trnaseq/modifications/__init__.py`
   - Implement `ModificationAnalyzer` class
   - Write docstrings

### Phase 2: CLI Interface (Week 2)

**Tasks:**
1. **Day 1-2:** CLI framework
   - Setup `click` framework in `trnaseq/cli/main.py`
   - Create subcommands structure
   - Add config file parsing (YAML)

2. **Day 3:** Modification CLI
   - Implement `trnaseq modifications` command
   - Add all options (--bam, --reference, --organism, etc.)
   - Test on sample data

3. **Day 4:** Fragment CLI
   - Implement `trnaseq fragments` commands
   - Add classify, quantify, differential subcommands

4. **Day 5:** Integration
   - Create `trnaseq run` for full pipeline
   - Add background execution
   - Test end-to-end

### Phase 3: Packaging & Distribution (Week 3)

**Tasks:**
1. **Day 1-2:** Packaging
   - Create `pyproject.toml` and `setup.py`
   - Add entry points for CLI
   - Create conda recipe

2. **Day 3:** Containerization
   - Write `Dockerfile`
   - Build and test Docker image
   - Create Singularity definition

3. **Day 4:** Nextflow modules
   - Create individual process modules
   - Test with example data

4. **Day 5:** Documentation
   - Setup Sphinx
   - Write quickstart guide
   - Create tutorial notebooks

---

## Testing Strategy

### Unit Tests

```python
# tests/unit/test_modifications.py
import pytest
from trnaseq.modifications import RTSignatureDetector

def test_rt_signature_detection():
    detector = RTSignatureDetector(min_coverage=50)
    signatures = detector.analyze_bam(
        'tests/data/sample.bam',
        'tests/data/ref.fa'
    )

    assert len(signatures) > 0
    assert 'position' in signatures.columns
    assert 'mismatch_rate' in signatures.columns

def test_modification_calling():
    from trnaseq.modifications import ModificationCaller

    caller = ModificationCaller(organism='human')
    # Test with synthetic data
    ...
```

### Integration Tests

```python
# tests/integration/test_full_pipeline.py
def test_modification_pipeline():
    """Test full modification analysis pipeline"""
    from trnaseq import ModificationAnalyzer

    analyzer = ModificationAnalyzer(
        bam_files=['tests/data/sample1.bam'],
        reference='tests/data/hg38-tRNAs.fa',
        organism='human'
    )

    signatures = analyzer.detect_rt_signatures()
    modifications = analyzer.call_modifications(signatures)

    assert len(modifications) > 0
    assert 'm1A' in modifications['modification_type'].values
```

### CLI Tests

```bash
# tests/cli/test_modifications.sh
#!/bin/bash

# Test modification detection
trnaseq modifications \
    --bam tests/data/sample1.bam \
    --reference tests/data/hg38-tRNAs.fa \
    --output tests/output/ \
    --signatures-only

# Check output exists
test -f tests/output/rt_signatures.csv || exit 1

# Test modification calling
trnaseq modifications call \
    --signatures tests/output/rt_signatures.csv \
    --organism human \
    --output tests/output/

test -f tests/output/modifications.csv || exit 1
```

---

## Migration Strategy

### Backwards Compatibility

**Goal:** Existing notebooks continue to work while new CLI is available

**Approach:**
1. Keep existing code in `src/` as-is
2. Create new package `trnaseq/` with refactored code
3. Add adapter layer for compatibility

```python
# trnaseq/compat.py - Backwards compatibility layer

from trnaseq.modifications import RTSignatureDetector

class TM_analysis:
    """
    Backwards-compatible wrapper for notebook code

    This allows existing notebooks using TM_analysis to work
    while internally using the new modular code.
    """
    def __init__(self, dir_dict, sample_df, tRNA_database, **kwargs):
        # Map old parameters to new API
        self._detector = RTSignatureDetector(**self._map_kwargs(kwargs))
        self.dir_dict = dir_dict
        self.sample_df = sample_df

    def find_muts(self, *args, **kwargs):
        """Old method name → new API"""
        return self._detector.detect_signatures(*args, **kwargs)

    # ... other compatibility methods
```

**Usage in existing notebooks:**
```python
# Old code (still works)
from src.transcript_mutations import TM_analysis

# New code (preferred)
from trnaseq import ModificationAnalyzer
```

---

## Documentation Structure

### User Documentation

1. **Installation Guide**
   - PyPI, conda, Docker, Singularity
   - Dependency requirements
   - HMS O2 setup

2. **Quickstart**
   - 5-minute tutorial
   - Example data download
   - Basic commands

3. **User Guide**
   - Preprocessing workflow
   - Alignment parameters
   - Modification analysis
   - Fragment analysis
   - Visualization options

4. **CLI Reference**
   - All commands and options
   - Examples for each command

5. **Notebook Tutorials**
   - Interactive Jupyter notebooks
   - Real-world examples
   - RVHMS04 and RVHMS09 workflows

### Developer Documentation

1. **API Reference**
   - Auto-generated from docstrings
   - Module descriptions
   - Class diagrams

2. **Contributing Guide**
   - Code style (Black, flake8)
   - Testing requirements
   - Pull request process

3. **Architecture**
   - Design decisions
   - Module dependencies
   - Extension points

---

## Dependency Management

### Core Dependencies

```toml
# pyproject.toml
[project]
name = "trnaseq-charge"
version = "1.0.0"
description = "Comprehensive tRNA-seq analysis toolkit"
authors = [{name = "Russel Vincent", email = "russel_vincent@hms.harvard.edu"}]
license = {text = "MIT"}
requires-python = ">=3.8"

dependencies = [
    "numpy>=1.24",
    "pandas>=2.0",
    "scipy>=1.10",
    "biopython>=1.81",
    "pysam>=0.21",
    "click>=8.0",
    "pyyaml>=6.0",
    "matplotlib>=3.7",
    "seaborn>=0.12",
]

[project.optional-dependencies]
plotting = [
    "logomaker",
    "adjustText",
]
ml = [
    "scikit-learn>=1.3",
    "statsmodels>=0.14",
]
nextflow = [
    "nextflow>=22.0",
]
dev = [
    "pytest>=7.0",
    "black>=23.0",
    "flake8>=6.0",
    "sphinx>=6.0",
]

[project.scripts]
trnaseq = "trnaseq.cli.main:cli"
```

---

## Continuous Integration

### GitHub Actions

```yaml
# .github/workflows/test.yml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.8", "3.9", "3.10", "3.11"]

    steps:
      - uses: actions/checkout@v3

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install dependencies
        run: |
          pip install -e ".[dev]"

      - name: Run tests
        run: |
          pytest tests/ --cov=trnaseq --cov-report=xml

      - name: Upload coverage
        uses: codecov/codecov-action@v3
```

---

## Next Steps

1. **Today:** Start extracting RT signature module
2. **This week:** Complete modification analysis module + CLI
3. **Next week:** Package for distribution
4. **Week 3:** Nextflow integration + documentation

**First PR:** Modification analysis module (Option B implementation)

---

**Last Updated:** 2026-01-23
