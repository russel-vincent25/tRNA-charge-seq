"""
RT Signature Analysis for tRNA Modifications

This module analyzes reverse transcription (RT) signatures in tRNA sequencing data
to identify positions with elevated mismatch rates, gaps, and RT stops - indicators
of post-transcriptional modifications.

Key concepts:
- PSCM (Position-Specific Count Matrix): Tracks nucleotide observations at each position
- RT signatures: Mismatches, deletions, or RT stops caused by modified nucleotides
- Modification sites: Positions with statistically elevated RT signature rates
"""

import numpy as np
import pandas as pd
import bz2
import warnings
from Bio import SeqIO, Align
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


class RTSignatureAnalyzer:
    """
    Analyze RT signatures from tRNA-seq alignment data.

    This class processes alignment data (BAM files or alignment statistics)
    to identify positions in tRNA transcripts that show elevated mismatch rates,
    gaps, or RT stops - all indicators of post-transcriptional modifications.

    Attributes:
        char_list: Nucleotide characters tracked (A, C, G, T, U, N, -)
        char_dict: Mapping of characters to indices
        min_coverage: Minimum coverage required to call a signature
        mismatch_threshold: Minimum mismatch rate to flag a position

    Example:
        >>> analyzer = RTSignatureAnalyzer(min_coverage=50)
        >>> signatures = analyzer.analyze_from_stats(
        ...     stats_csv='sample1_stats.csv.bz2',
        ...     reference_fasta='hg38-tRNAs.fa'
        ... )
        >>> print(signatures.head())
    """

    def __init__(
        self,
        min_coverage: int = 50,
        mismatch_threshold: float = 0.10,
        rt_stop_threshold: float = 20.0,
        verbose: bool = True
    ):
        """
        Initialize RT signature analyzer.

        Args:
            min_coverage: Minimum read coverage to analyze a position (default: 50)
            mismatch_threshold: Minimum mismatch rate to flag as signature (default: 0.10 = 10%)
            rt_stop_threshold: Minimum RT stop percentage to flag (default: 20.0%)
            verbose: Print progress messages (default: True)
        """
        self.min_coverage = min_coverage
        self.mismatch_threshold = mismatch_threshold
        self.rt_stop_threshold = rt_stop_threshold
        self.verbose = verbose

        # Nucleotide character tracking
        self.char_str = 'ACGTUN-'
        self.char_list = [c for c in self.char_str]
        self.char_dict = {c: i for i, c in enumerate(self.char_str)}

        # Storage for analysis results
        self.reference_sequences = {}
        self.pscm_data = {}  # Position-Specific Count Matrices

    def load_reference(self, reference_fasta: Union[str, Path]) -> Dict[str, Dict]:
        """
        Load reference tRNA sequences from FASTA file.

        Args:
            reference_fasta: Path to FASTA file with tRNA sequences

        Returns:
            Dictionary mapping tRNA names to sequence info (seq, seq_len)
        """
        ref_dict = {}
        for record in SeqIO.parse(str(reference_fasta), "fasta"):
            ref_dict[record.id] = {
                'seq': str(record.seq),
                'seq_len': len(record.seq)
            }

        if self.verbose:
            print(f"Loaded {len(ref_dict)} reference tRNA sequences")

        self.reference_sequences = ref_dict
        return ref_dict

    def initialize_pscm(self, trna_name: str) -> np.ndarray:
        """
        Initialize Position-Specific Count Matrix for a tRNA.

        Args:
            trna_name: Name of the tRNA transcript

        Returns:
            Zero-initialized PSCM matrix (seq_len × 7 nucleotides)
        """
        if trna_name not in self.reference_sequences:
            raise ValueError(f"tRNA {trna_name} not found in reference sequences")

        seq_len = self.reference_sequences[trna_name]['seq_len']
        return np.zeros((seq_len, len(self.char_list)))

    def build_pscm_from_alignment(
        self,
        read_seq: str,
        reference_seq: str,
        weight: float = 1.0,
        aligner: Optional[Align.PairwiseAligner] = None
    ) -> np.ndarray:
        """
        Build PSCM for a single read-reference alignment.

        .. deprecated::
            Use :class:`~trnaseq.modifications.positional.PositionalExtractor`
            for efficient batch PSCM extraction from SWalign JSON files.

        Performs pairwise alignment between read and reference, then extracts
        the observed nucleotide at each position (including mismatches and gaps).

        Args:
            read_seq: Read sequence
            reference_seq: Reference tRNA sequence
            weight: Weight for this observation (default: 1.0, for UMI counts use fractional)
            aligner: Pairwise aligner object (creates default if None)

        Returns:
            PSCM matrix for this alignment (seq_len × 7)
        """
        # Create aligner if not provided
        if aligner is None:
            aligner = Align.PairwiseAligner()
            aligner.mode = 'local'
            aligner.match_score = 1
            aligner.mismatch_score = -2
            aligner.open_gap_score = -3
            aligner.extend_gap_score = -2

        # Perform alignment
        alignments = aligner.align(reference_seq, read_seq)

        # Take first (best) alignment
        if len(alignments) == 0:
            warnings.warn(f"No alignment found for read")
            return np.zeros((len(reference_seq), len(self.char_list)))

        alignment = alignments[0]
        t_cor, q_cor = alignment.aligned

        # Initialize count matrix
        count_mat = np.zeros((len(reference_seq), len(self.char_list)))

        # Find gaps in reference (deletions in read)
        for i in range(1, len(t_cor)):
            for j in range(t_cor[i][0] - t_cor[i-1][1]):
                gap_pos = t_cor[i-1][1] + j
                count_mat[gap_pos, self.char_dict['-']] = weight

        # Find matches and mismatches
        for tran, qran in zip(t_cor, q_cor):
            for ti, qi in zip(range(*tran), range(*qran)):
                observed_nt = read_seq[qi]
                count_mat[ti, self.char_dict[observed_nt]] = weight

        return count_mat

    def calculate_mismatch_rates(
        self,
        pscm_df: pd.DataFrame,
        reference_seq: str,
        min_coverage: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Calculate per-position mismatch rates from PSCM.

        Mismatch rate = (total_coverage - correct_nt_count) / total_coverage

        Args:
            pscm_df: DataFrame with PSCM data (positions × nucleotides)
            reference_seq: Reference sequence to determine correct nucleotide
            min_coverage: Minimum coverage to include (uses instance default if None)

        Returns:
            DataFrame with columns: position, coverage, mismatches, mismatch_rate, correct_nt
        """
        if min_coverage is None:
            min_coverage = self.min_coverage

        # Calculate total coverage per position
        coverage = pscm_df.sum(axis=1).values

        # Determine correct nucleotide at each position
        results = []
        for pos in range(len(reference_seq)):
            if coverage[pos] < min_coverage:
                continue

            correct_nt = reference_seq[pos]
            if correct_nt not in self.char_dict:
                continue

            correct_count = pscm_df.loc[pos, correct_nt]
            mismatch_count = coverage[pos] - correct_count
            mismatch_rate = mismatch_count / coverage[pos]

            results.append({
                'position': pos + 1,  # 1-based positioning
                'coverage': coverage[pos],
                'correct_nt': correct_nt,
                'correct_count': correct_count,
                'mismatch_count': mismatch_count,
                'mismatch_rate': mismatch_rate
            })

        return pd.DataFrame(results)

    def calculate_gap_rates(
        self,
        pscm_df: pd.DataFrame,
        min_coverage: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Calculate per-position gap (deletion) rates from PSCM.

        Gap rate = gap_count / total_coverage

        Args:
            pscm_df: DataFrame with PSCM data (positions × nucleotides)
            min_coverage: Minimum coverage to include (uses instance default if None)

        Returns:
            DataFrame with columns: position, coverage, gap_count, gap_rate
        """
        if min_coverage is None:
            min_coverage = self.min_coverage

        # Calculate total coverage per position
        coverage = pscm_df.sum(axis=1).values
        gap_counts = pscm_df['-'].values

        results = []
        for pos in range(len(coverage)):
            if coverage[pos] < min_coverage:
                continue

            gap_rate = gap_counts[pos] / coverage[pos]

            results.append({
                'position': pos + 1,  # 1-based
                'coverage': coverage[pos],
                'gap_count': gap_counts[pos],
                'gap_rate': gap_rate
            })

        return pd.DataFrame(results)

    def calculate_rt_stops(
        self,
        pscm_df: pd.DataFrame,
        rt_stop_counts: Optional[np.ndarray] = None,
    ) -> pd.DataFrame:
        """
        Calculate RT stop frequencies.

        RT stops indicate positions where reverse transcriptase prematurely
        terminates, often due to modified nucleotides (e.g., m1A).

        When ``rt_stop_counts`` is provided (from PositionalExtractor), the
        actual per-position RT stop counts are used directly. Otherwise, the
        heuristic coverage-difference formula is applied:

            RTstop(i) = 100 * (cov(i+1) - cov(i)) / cov(i+1)

        Args:
            pscm_df: DataFrame with PSCM data (positions x nucleotides)
            rt_stop_counts: Optional array of direct RT stop counts per
                position (from PositionalExtractor PSCM column 7).

        Returns:
            DataFrame with columns: position, coverage, coverage_next,
            rt_stop_pct, rt_stop_count (when available)
        """
        # Calculate total coverage per position
        coverage = pscm_df.sum(axis=1).values

        if rt_stop_counts is not None:
            # Use actual RT stop counts from PositionalExtractor
            rt_stop_counts = np.asarray(rt_stop_counts, dtype=np.float64)
            rt_stops = np.divide(
                rt_stop_counts * 100,
                coverage,
                out=np.zeros_like(coverage, dtype=np.float64),
                where=coverage != 0,
            )
            coverage_next = np.roll(coverage, -1)
            coverage_next[-1] = coverage[-1]

            results = []
            for pos in range(len(coverage)):
                results.append({
                    'position': pos + 1,
                    'coverage': coverage[pos],
                    'coverage_next': coverage_next[pos],
                    'rt_stop_pct': rt_stops[pos],
                    'rt_stop_count': int(rt_stop_counts[pos]),
                })
        else:
            # Heuristic: coverage-difference formula (Wang et al. 2021)
            coverage_next = np.roll(coverage, -1)
            coverage_next[-1] = coverage[-1]

            rt_stops = 100 * np.divide(
                (coverage_next - coverage),
                coverage_next,
                out=np.zeros_like(coverage, dtype=np.float64),
                where=coverage_next != 0
            )

            results = []
            for pos in range(len(coverage)):
                results.append({
                    'position': pos + 1,
                    'coverage': coverage[pos],
                    'coverage_next': coverage_next[pos],
                    'rt_stop_pct': rt_stops[pos],
                })

        return pd.DataFrame(results)

    def identify_signature_positions(
        self,
        mismatch_df: pd.DataFrame,
        gap_df: pd.DataFrame,
        rt_stop_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Identify positions with significant RT signatures.

        Combines mismatch, gap, and RT stop data to flag positions that
        show elevated rates indicative of modifications.

        Args:
            mismatch_df: Mismatch rate DataFrame
            gap_df: Gap rate DataFrame
            rt_stop_df: RT stop DataFrame

        Returns:
            DataFrame with flagged signature positions and evidence
        """
        # Guard against empty DataFrames (all positions below min_coverage)
        _gap = (
            gap_df[['position', 'gap_rate']]
            if not gap_df.empty and 'position' in gap_df.columns
            else pd.DataFrame(columns=['position', 'gap_rate'])
        )
        _rt = (
            rt_stop_df[['position', 'rt_stop_pct']]
            if not rt_stop_df.empty and 'position' in rt_stop_df.columns
            else pd.DataFrame(columns=['position', 'rt_stop_pct'])
        )
        if mismatch_df.empty or 'position' not in mismatch_df.columns:
            mismatch_df = pd.DataFrame(columns=[
                'position', 'coverage', 'correct_nt', 'correct_count',
                'mismatch_count', 'mismatch_rate',
            ])

        # Merge all data
        merged = mismatch_df.merge(
            _gap, on='position', how='outer'
        ).merge(
            _rt, on='position', how='outer'
        )

        # Fill NaN values
        merged = merged.fillna(0)

        # If completely empty, return with expected columns
        if merged.empty:
            for col in ['mismatch_rate', 'rt_stop_pct', 'gap_rate',
                        'has_mismatch_signature', 'has_rt_stop',
                        'has_gap_signature', 'has_signature']:
                merged[col] = pd.Series(dtype='float64')
            return merged

        # Flag positions based on thresholds
        merged['has_mismatch_signature'] = merged['mismatch_rate'] >= self.mismatch_threshold
        merged['has_rt_stop'] = merged['rt_stop_pct'] >= self.rt_stop_threshold
        merged['has_gap_signature'] = merged['gap_rate'] >= 0.05  # 5% gap threshold

        # Overall signature flag
        merged['has_signature'] = (
            merged['has_mismatch_signature'] |
            merged['has_rt_stop'] |
            merged['has_gap_signature']
        )

        # Sort by position
        merged = merged.sort_values('position')

        return merged

    def extract_mutation_patterns(
        self,
        pscm_df: pd.DataFrame,
        position: int,
        reference_seq: str
    ) -> Dict[str, float]:
        """
        Extract detailed mutation patterns at a specific position.

        Args:
            pscm_df: DataFrame with PSCM data
            position: Position to analyze (1-based)
            reference_seq: Reference sequence

        Returns:
            Dictionary of mutation types and counts (e.g., {'A->G': 45, 'A->C': 3})
        """
        pos_idx = position - 1  # Convert to 0-based

        if pos_idx < 0 or pos_idx >= len(reference_seq):
            raise ValueError(f"Position {position} out of range")

        ref_nt = reference_seq[pos_idx]
        patterns = {}

        # Get counts for each observed nucleotide
        for nt in self.char_list:
            count = pscm_df.loc[pos_idx, nt]
            if count > 0:
                if nt == '-':
                    patterns['deletion'] = count
                elif nt == ref_nt:
                    patterns['match'] = count
                else:
                    patterns[f'{ref_nt}->{nt}'] = count

        return patterns

    def analyze_trna(
        self,
        trna_name: str,
        pscm_df: pd.DataFrame
    ) -> Dict[str, pd.DataFrame]:
        """
        Complete RT signature analysis for a single tRNA.

        Args:
            trna_name: Name of the tRNA
            pscm_df: Position-Specific Count Matrix as DataFrame

        Returns:
            Dictionary with analysis results:
                - 'mismatch': Mismatch rate DataFrame
                - 'gaps': Gap rate DataFrame
                - 'rt_stops': RT stop DataFrame
                - 'signatures': Combined signature positions
        """
        if trna_name not in self.reference_sequences:
            raise ValueError(f"tRNA {trna_name} not found in reference")

        ref_seq = self.reference_sequences[trna_name]['seq']

        # Calculate all metrics
        mismatch_df = self.calculate_mismatch_rates(pscm_df, ref_seq)
        gap_df = self.calculate_gap_rates(pscm_df)
        rt_stop_df = self.calculate_rt_stops(pscm_df)

        # Identify signature positions
        signatures_df = self.identify_signature_positions(
            mismatch_df, gap_df, rt_stop_df
        )

        return {
            'mismatch': mismatch_df,
            'gaps': gap_df,
            'rt_stops': rt_stop_df,
            'signatures': signatures_df
        }

    def process_sample_from_pipeline(
        self,
        stats_csv: Union[str, Path],
        umi_trimmed_fastq: Union[str, Path],
        species: str,
        use_umi_count: bool = True,
        unique_anno: bool = True,
        max_5p_non_temp: int = 10
    ) -> Dict[str, Dict]:
        """
        Process a sample using the existing pipeline's output format.

        .. deprecated::
            Use :class:`~trnaseq.modifications.positional.PositionalExtractor`
            which streams SWalign JSON files directly for much better performance.

        This method is compatible with the existing tRNA-charge-seq pipeline,
        reading stats CSV and UMI-trimmed FASTQ files to build PSCMs.

        Args:
            stats_csv: Path to stats CSV file (can be .bz2 compressed)
            umi_trimmed_fastq: Path to UMI-trimmed FASTQ file (can be .bz2 compressed)
            species: Species name (must match reference key)
            use_umi_count: Use UMI counts instead of read counts
            unique_anno: Only use reads with unique annotation
            max_5p_non_temp: Maximum allowed 5' non-template length

        Returns:
            Dictionary mapping tRNA names to their PSCM DataFrames
        """
        count_col = 'UMIcount' if use_umi_count else 'count'

        # Read stats CSV
        if str(stats_csv).endswith('.bz2'):
            with bz2.open(stats_csv, 'rt', encoding='utf-8') as fh:
                stats_df = pd.read_csv(fh, keep_default_na=False)
        else:
            stats_df = pd.read_csv(stats_csv, keep_default_na=False)

        # Filter stats based on quality criteria
        mask = (
            stats_df['3p_cover'] &
            (stats_df['3p_non-temp'] == '') &
            (stats_df['5p_non-temp'].apply(len) <= max_5p_non_temp) &
            ((stats_df['align_3p_nt'] == 'A') | (stats_df['align_3p_nt'] == 'C'))
        )
        stats_df = stats_df[mask]

        # Build readID to annotation mapping
        id_to_anno = {}
        for _, row in stats_df.iterrows():
            id_to_anno[row['readID']] = {
                'count': row[count_col],
                'anno': row['tRNA_annotation'].split('@')
            }

        # Deduplicate sequences from FASTQ
        if str(umi_trimmed_fastq).endswith('.bz2'):
            fh = bz2.open(umi_trimmed_fastq, 'rt')
        else:
            fh = open(umi_trimmed_fastq, 'r')

        dedup_seq_count = {}
        for record in SeqIO.parse(fh, 'fastq'):
            seq = str(record.seq)
            if seq in dedup_seq_count:
                dedup_seq_count[seq]['count'] += 1
            else:
                dedup_seq_count[seq] = {
                    'count': 1,
                    'id': record.id
                }
        fh.close()

        # Initialize PSCMs for all tRNAs
        pscm_results = {}
        for trna_name in self.reference_sequences:
            seq_len = self.reference_sequences[trna_name]['seq_len']
            pscm_results[trna_name] = np.zeros((seq_len, len(self.char_list)))

        # Create aligner
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 1
        aligner.mismatch_score = -2
        aligner.open_gap_score = -3
        aligner.extend_gap_score = -2

        # Process each unique sequence
        for seq, seq_info in dedup_seq_count.items():
            read_id = seq_info['id']

            if read_id not in id_to_anno:
                continue

            anno_list = id_to_anno[read_id]['anno']
            seq_count = id_to_anno[read_id]['count']

            # Skip if multiple annotations and unique requested
            if unique_anno and len(anno_list) > 1:
                continue

            # Generate alignments for all annotations
            alignments_anno = []
            for anno in anno_list:
                if anno not in self.reference_sequences:
                    continue
                target = self.reference_sequences[anno]['seq']
                alignments = aligner.align(target, seq)
                alignments_anno.append((anno, alignments))

            # Weight by number of equivalent alignments
            total_alignments = sum(len(alns) for _, alns in alignments_anno)
            if total_alignments == 0:
                continue

            weight = 1.0 / total_alignments * seq_count

            # Process each alignment
            for anno, alignments in alignments_anno:
                target = self.reference_sequences[anno]['seq']
                for alignment in alignments:
                    # Build PSCM for this alignment
                    count_mat = self.build_pscm_from_alignment(
                        seq, target, weight=weight, aligner=None
                    )
                    pscm_results[anno] += count_mat

        # Convert to DataFrames
        pscm_dfs = {}
        for trna_name, pscm_array in pscm_results.items():
            if pscm_array.sum() > 0:  # Only include tRNAs with data
                pscm_dfs[trna_name] = pd.DataFrame(pscm_array, columns=self.char_list)

        if self.verbose:
            print(f"Processed {len(dedup_seq_count)} unique sequences")
            print(f"Generated PSCMs for {len(pscm_dfs)} tRNAs")

        self.pscm_data = pscm_dfs
        return pscm_dfs

    def load_pscm_from_positional(
        self,
        pscm_dict: Dict[str, 'np.ndarray'],
    ) -> Dict[str, pd.DataFrame]:
        """Load pre-computed PSCMs from PositionalExtractor.

        Converts PositionalExtractor's ndarray format (ref_len x 8) into the
        DataFrame format expected by this class (ref_len x 7, columns ACGTUN-).

        Args:
            pscm_dict: Output of ``PositionalExtractor.extract_sample()``.
                Keys are tRNA names, values are ndarrays with columns
                [A, C, G, T, N, gap, coverage, rt_stop].

        Returns:
            Dictionary mapping tRNA names to PSCM DataFrames compatible
            with all analysis methods in this class.
        """
        pscm_dfs = {}
        for trna_name, mat in pscm_dict.items():
            # Map PositionalExtractor columns to RTSignatureAnalyzer columns
            # Extractor: [A=0, C=1, G=2, T=3, N=4, gap=5, coverage=6, rt_stop=7]
            # Analyzer:  [A, C, G, T, U, N, -]
            n_pos = mat.shape[0]
            df_data = {
                'A': mat[:, 0],
                'C': mat[:, 1],
                'G': mat[:, 2],
                'T': mat[:, 3],
                'U': np.zeros(n_pos),  # U not tracked separately
                'N': mat[:, 4],
                '-': mat[:, 5],
            }
            pscm_dfs[trna_name] = pd.DataFrame(df_data)

        if self.verbose:
            print(f"Loaded PSCMs for {len(pscm_dfs)} tRNAs from PositionalExtractor")

        self.pscm_data = pscm_dfs
        # Also store raw matrices for rt_stop counts
        self._positional_raw = pscm_dict
        return pscm_dfs

    def analyze_trna_with_actual_stops(
        self,
        trna_name: str,
        pscm_df: pd.DataFrame,
        rt_stop_counts: Optional[np.ndarray] = None,
    ) -> Dict[str, pd.DataFrame]:
        """Analyze a single tRNA using actual RT stop counts when available.

        Like :meth:`analyze_trna`, but passes RT stop counts directly to
        :meth:`calculate_rt_stops` for more accurate profiling.

        Args:
            trna_name: Name of the tRNA.
            pscm_df: Position-Specific Count Matrix as DataFrame.
            rt_stop_counts: Per-position RT stop counts from PositionalExtractor.

        Returns:
            Dictionary with analysis results (mismatch, gaps, rt_stops, signatures).
        """
        if trna_name not in self.reference_sequences:
            raise ValueError(f"tRNA {trna_name} not found in reference")

        ref_seq = self.reference_sequences[trna_name]['seq']

        mismatch_df = self.calculate_mismatch_rates(pscm_df, ref_seq)
        gap_df = self.calculate_gap_rates(pscm_df)
        rt_stop_df = self.calculate_rt_stops(pscm_df, rt_stop_counts=rt_stop_counts)

        signatures_df = self.identify_signature_positions(
            mismatch_df, gap_df, rt_stop_df
        )

        return {
            'mismatch': mismatch_df,
            'gaps': gap_df,
            'rt_stops': rt_stop_df,
            'signatures': signatures_df,
        }

    def analyze_all_trnas(
        self,
        pscm_dict: Optional[Dict[str, pd.DataFrame]] = None
    ) -> Dict[str, Dict[str, pd.DataFrame]]:
        """
        Analyze RT signatures for all tRNAs.

        Args:
            pscm_dict: Dictionary of tRNA name -> PSCM DataFrame.
                      If None, uses self.pscm_data

        Returns:
            Dictionary mapping tRNA names to their analysis results
        """
        if pscm_dict is None:
            pscm_dict = self.pscm_data

        if not pscm_dict:
            raise ValueError("No PSCM data available. Run process_sample_from_pipeline first.")

        results = {}
        for trna_name, pscm_df in pscm_dict.items():
            try:
                results[trna_name] = self.analyze_trna(trna_name, pscm_df)
            except Exception as e:
                if self.verbose:
                    print(f"Warning: Failed to analyze {trna_name}: {e}")
                continue

        return results


# Utility function for quick analysis
def analyze_rt_signatures(
    reference_fasta: Union[str, Path],
    min_coverage: int = 50,
    mismatch_threshold: float = 0.10,
    verbose: bool = True
) -> RTSignatureAnalyzer:
    """
    Quick initialization of RT signature analyzer with reference loading.

    Args:
        reference_fasta: Path to reference tRNA FASTA
        min_coverage: Minimum coverage threshold
        mismatch_threshold: Mismatch rate threshold
        verbose: Print progress

    Returns:
        Initialized RTSignatureAnalyzer with loaded references
    """
    analyzer = RTSignatureAnalyzer(
        min_coverage=min_coverage,
        mismatch_threshold=mismatch_threshold,
        verbose=verbose
    )
    analyzer.load_reference(reference_fasta)
    return analyzer
