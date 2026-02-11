# tRNA Modification Analysis Module

Standalone module for detecting and analyzing tRNA modifications based on RT (reverse transcription) signature analysis.

## Overview

This module extracts and modernizes the RT signature analysis code from the original `transcript_mutations.py` into a clean, reusable package. It identifies positions in tRNA transcripts that show elevated mismatch rates, gaps, or RT stops - all indicators of post-transcriptional modifications.

## Installation

```python
# From the tRNA-charge-seq directory
from trnaseq.modifications import RTSignatureAnalyzer
```

## Quick Start

```python
from trnaseq.modifications import RTSignatureAnalyzer

# Initialize analyzer
analyzer = RTSignatureAnalyzer(
    min_coverage=50,
    mismatch_threshold=0.10,
    rt_stop_threshold=20.0
)

# Load reference tRNA sequences
analyzer.load_reference('path/to/hg38-tRNAs.fa')

# Process a sample from the existing pipeline
pscm_dict = analyzer.process_sample_from_pipeline(
    stats_csv='data/stats_collection/sample1_stats.csv.bz2',
    umi_trimmed_fastq='data/UMI_trimmed/sample1_UMI-trimmed.fastq.bz2',
    species='human',
    use_umi_count=True
)

# Analyze all tRNAs
results = analyzer.analyze_all_trnas()

# Inspect results
for trna_name, trna_results in results.items():
    signatures = trna_results['signatures']
    print(f"\n{trna_name}:")
    print(signatures[signatures['has_signature']])
```

## Key Features

### RT Signature Detection

The module identifies three types of RT signatures:

1. **Mismatch signatures**: Elevated mismatch rates at specific positions
   - Example: m1A at position 58 causes A→G, A→C mismatches

2. **Gap signatures**: Deletions in reads at modification sites
   - Example: m3C can cause deletion signatures

3. **RT stop signatures**: Premature reverse transcriptase termination
   - Example: m1A causes RT stops before position 58

### Position-Specific Count Matrix (PSCM)

Tracks nucleotide observations (A, C, G, T, U, N, -) at each position:
- Built from read alignments to reference sequences
- Weighted by UMI counts or read counts
- Accounts for alignment ambiguity

### Compatibility

- **Backward compatible** with existing tRNA-charge-seq pipeline output
- Reads stats CSV files and UMI-trimmed FASTQ files
- Can also work with BAM files (future enhancement)

## Module Components

### `rt_signatures.py`

Core module containing:

- `RTSignatureAnalyzer`: Main analysis class
- Methods for:
  - Loading reference sequences
  - Building PSCMs from alignments
  - Calculating mismatch rates
  - Calculating gap rates
  - Calculating RT stops
  - Identifying signature positions
  - Processing pipeline output files

### `example_usage.py`

Example scripts demonstrating:
- Basic RT signature analysis
- Single tRNA analysis
- Exporting results to CSV
- Integration with existing pipeline

## API Reference

### RTSignatureAnalyzer

```python
analyzer = RTSignatureAnalyzer(
    min_coverage=50,          # Minimum coverage to analyze position
    mismatch_threshold=0.10,  # Minimum mismatch rate for signature (10%)
    rt_stop_threshold=20.0,   # Minimum RT stop percentage for signature
    verbose=True              # Print progress messages
)
```

### Key Methods

#### `load_reference(reference_fasta)`
Load reference tRNA sequences from FASTA file.

#### `process_sample_from_pipeline(stats_csv, umi_trimmed_fastq, species, ...)`
Process a sample using existing pipeline output format.

#### `analyze_trna(trna_name, pscm_df)`
Complete RT signature analysis for a single tRNA.

Returns dictionary with:
- `'mismatch'`: Mismatch rate DataFrame
- `'gaps'`: Gap rate DataFrame
- `'rt_stops'`: RT stop DataFrame
- `'signatures'`: Combined signature positions

#### `analyze_all_trnas(pscm_dict=None)`
Analyze RT signatures for all tRNAs.

#### `extract_mutation_patterns(pscm_df, position, reference_seq)`
Get detailed mutation patterns at a specific position.

## Output Format

### Signature DataFrame

Each tRNA's signature analysis returns a DataFrame with columns:

- `position`: 1-based position in tRNA
- `coverage`: Read coverage at position
- `correct_nt`: Reference nucleotide
- `mismatch_count`: Number of mismatches
- `mismatch_rate`: Fraction of reads with mismatches
- `gap_rate`: Fraction of reads with gaps
- `rt_stop_pct`: Percentage RT stop signal
- `has_mismatch_signature`: Boolean flag (>10% mismatch)
- `has_rt_stop`: Boolean flag (>20% RT stop)
- `has_gap_signature`: Boolean flag (>5% gap)
- `has_signature`: Overall signature flag (any of the above)

## Examples

### Export High-Confidence Modification Sites

```python
import pandas as pd

# Run analysis
analyzer = RTSignatureAnalyzer()
analyzer.load_reference('hg38-tRNAs.fa')
pscm_dict = analyzer.process_sample_from_pipeline(...)
results = analyzer.analyze_all_trnas()

# Combine all signatures
all_sigs = []
for trna_name, res in results.items():
    df = res['signatures'].copy()
    df['trna_name'] = trna_name
    all_sigs.append(df)

combined = pd.concat(all_sigs, ignore_index=True)

# Filter high-confidence sites
high_conf = combined[
    (combined['mismatch_rate'] > 0.15) |  # >15% mismatch
    (combined['rt_stop_pct'] > 30)         # >30% RT stop
]

high_conf.to_csv('high_confidence_modifications.csv', index=False)
```

### Analyze Specific Modification Site

```python
# Known m1A site at position 58 in tRNA-Thr
trna_name = 'tRNA-Thr-AGT-1-1'
results = analyzer.analyze_trna(trna_name, pscm_dict[trna_name])

# Check position 58
sig_df = results['signatures']
pos58 = sig_df[sig_df['position'] == 58]

print(f"Position 58 mismatch rate: {pos58['mismatch_rate'].values[0]:.2%}")
print(f"Position 58 RT stop: {pos58['rt_stop_pct'].values[0]:.1f}%")

# Get mutation patterns
patterns = analyzer.extract_mutation_patterns(
    pscm_dict[trna_name],
    position=58,
    reference_seq=analyzer.reference_sequences[trna_name]['seq']
)
print("Mutation patterns:", patterns)
```

## Modification Calling (NEW!)

The module now includes automated modification calling based on RT signature patterns.

### Quick Start with Modification Calling

```python
from trnaseq.modifications import RTSignatureAnalyzer, ModificationCaller

# Step 1: RT signature analysis
analyzer = RTSignatureAnalyzer(min_coverage=50)
analyzer.load_reference('hg38-tRNAs.fa')
pscm_dict = analyzer.process_sample_from_pipeline(
    stats_csv='sample_stats.csv.bz2',
    umi_trimmed_fastq='sample_UMI-trimmed.fastq.bz2',
    species='human'
)
rt_results = analyzer.analyze_all_trnas()

# Step 2: Call modifications
caller = ModificationCaller(
    organism='human',
    min_confidence=0.5,
    statistical_test=True
)

ref_seqs = {n: i['seq'] for n, i in analyzer.reference_sequences.items()}
modifications = caller.call_modifications_for_all_trnas(
    rt_results,
    pscm_dict,
    ref_seqs
)

# Display results
print(modifications[['trna_name', 'position', 'modification', 'confidence']])
```

### Supported Modifications

| Modification | Full Name | Typical Positions | RT Signature |
|--------------|-----------|-------------------|--------------|
| **m1A** | 1-methyladenosine | 58, 14, 9 | A→any mismatches + RT stops |
| **m3C** | 3-methylcytosine | 32 | C→T mismatches |
| **Ψ** | pseudouridine | 27, 28, 31, 32, 39, 40, 55 | U→C mismatches |
| **m7G** | 7-methylguanosine | 46 | G→A mismatches |
| **m5C** | 5-methylcytosine | 48, 49, 34, 40 | C→T mismatches (subtle) |
| **i6A** | N6-isopentenyladenosine | 37 | A→G mismatches |
| **m2G** | N2-methylguanosine | 10, 26 | G→A mismatches |
| **m22G** | N2,N2-dimethylguanosine | 26 | G→A mismatches |

### Confidence Scoring

Confidence scores (0-1) are calculated based on:

1. **Mismatch rate** - Higher rates increase confidence
2. **RT stop percentage** - Required for m1A, boosts confidence when present
3. **Pattern specificity** - Matches expected mutation pattern
4. **Position priors** - Higher confidence at known positions
5. **Coverage** - More reads = more reliable
6. **Statistical significance** - Binomial test for elevated mismatch rate

### Modification Caller API

```python
caller = ModificationCaller(
    organism='human',           # Organism for position priors
    min_confidence=0.5,         # Minimum confidence to report
    use_position_priors=True,   # Boost confidence at known positions
    statistical_test=True,      # Perform binomial test
    alpha=0.01                  # Significance level
)
```

#### Key Methods

**`call_modifications_for_trna(trna_name, signatures_df, pscm_df, reference_seq)`**
Call modifications for a single tRNA.

**`call_modifications_for_all_trnas(rt_signature_results, pscm_dict, reference_sequences)`**
Call modifications for all tRNAs at once.

**`filter_by_confidence(calls_df, min_confidence)`**
Filter calls by confidence threshold.

**`summarize_modifications(calls_df)`**
Summarize calls by modification type.

### Output Format

Modification calls DataFrame includes:

- `trna_name`: tRNA name
- `position`: Position in tRNA (1-based)
- `modification`: Modification type (e.g., 'm1A')
- `full_name`: Full chemical name
- `confidence`: Confidence score (0-1)
- `mismatch_rate`: Observed mismatch rate
- `rt_stop_pct`: RT stop percentage
- `gap_rate`: Gap rate
- `coverage`: Read coverage
- `pattern_fraction`: Fraction matching expected pattern
- `pvalue`: P-value from binomial test (if enabled)
- `in_typical_position`: Boolean flag for known positions

### Examples

See `example_modification_calling.py` for:
- Complete workflow from RT signatures to modification calls
- Single tRNA detailed analysis
- Multi-sample comparison
- Validation against known sites

## Next Steps

Remaining enhancements:

1. ✅ ~~**Modification Caller** (Task #2)~~ - **COMPLETE**
2. **Visualization Module** (Task #3): Publication-quality plots of RT signatures and modifications
3. **MODOMICS Integration** (Task #4): Cross-reference with known modification databases
4. **CLI Interface** (Task #5): Command-line access for HPC workflows

## References

- Wang et al. 2021 - RT stop methodology
- Original `transcript_mutations.py` - Source implementation
- MODOMICS database - tRNA modification reference

## Changes from Original Code

### Improvements

✅ **Modular design**: Clean separation from notebook workflow
✅ **Type hints**: Full type annotations for better IDE support
✅ **Docstrings**: Comprehensive documentation for all methods
✅ **Flexible input**: Can work with different data formats
✅ **Pandas output**: Easy-to-analyze DataFrame results
✅ **Error handling**: Better error messages and warnings

### Maintained Features

✅ **PSCM calculation**: Identical to original implementation
✅ **Alignment logic**: Same Biopython alignment approach
✅ **Weight handling**: Proper weighting for multi-mapping reads
✅ **UMI support**: Full UMI count integration
✅ **Quality filtering**: Same read filtering criteria

## Authors

- **Russel Vincent** - Module extraction and modernization
- **Original implementation** - tRNA-charge-seq `transcript_mutations.py`

## License

MIT (same as tRNA-charge-seq repository)
