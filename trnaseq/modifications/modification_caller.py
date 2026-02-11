"""
Modification Caller for tRNA RT Signatures

This module annotates RT signatures with specific tRNA modification types based on
known RT signature patterns from the literature.

Key modifications detected:
- m1A (1-methyladenosine): Strong RT stops, A->any mismatches at positions 58, 14
- m3C (3-methylcytosine): C->T mismatches at positions 32
- Ψ (pseudouridine): U->C signature at multiple positions
- m7G (7-methylguanosine): G->A signature at position 46
- m5C (5-methylcytosine): Subtle C->T signature at positions 48, 49
- i6A (N6-isopentenyladenosine): A->G at position 37

References:
- Carlile et al. 2014 (Nature) - Pseudouridine detection
- Schwartz et al. 2014 (Cell) - m1A, m3C detection
- Hauenschild et al. 2015 (NAR) - Modification signatures
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from scipy.stats import binomtest


@dataclass
class ModificationProfile:
    """
    Profile for a specific tRNA modification type.

    Attributes:
        name: Modification name (e.g., 'm1A')
        full_name: Full chemical name
        typical_positions: Common positions where this modification occurs
        signature_type: Type of RT signature ('mismatch', 'rt_stop', 'gap', 'combined')
        mismatch_pattern: Expected mismatch pattern (e.g., 'A->G', 'A->any')
        min_rate: Minimum rate to call this modification
        rt_stop_required: Whether RT stops are required for calling
        min_rt_stop_pct: Minimum RT stop percentage if required
    """
    name: str
    full_name: str
    typical_positions: List[int]
    signature_type: str
    mismatch_pattern: Optional[str] = None
    min_rate: float = 0.10
    rt_stop_required: bool = False
    min_rt_stop_pct: float = 15.0
    confidence_weight: float = 1.0


# Known modification profiles based on literature
MODIFICATION_PROFILES = {
    'm1A': ModificationProfile(
        name='m1A',
        full_name='1-methyladenosine',
        typical_positions=[58, 14, 9],
        signature_type='combined',
        mismatch_pattern='A->any',
        min_rate=0.12,
        rt_stop_required=True,
        min_rt_stop_pct=20.0,
        confidence_weight=1.5  # High confidence when both mismatch + RT stop
    ),

    'm3C': ModificationProfile(
        name='m3C',
        full_name='3-methylcytosine',
        typical_positions=[32],
        signature_type='mismatch',
        mismatch_pattern='C->T',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.2
    ),

    'pseudouridine': ModificationProfile(
        name='Ψ',
        full_name='pseudouridine',
        typical_positions=[27, 28, 31, 32, 39, 40, 55, 13, 38],
        signature_type='mismatch',
        mismatch_pattern='U->C',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm7G': ModificationProfile(
        name='m7G',
        full_name='7-methylguanosine',
        typical_positions=[46],
        signature_type='mismatch',
        mismatch_pattern='G->A',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.1
    ),

    'm5C': ModificationProfile(
        name='m5C',
        full_name='5-methylcytosine',
        typical_positions=[48, 49, 34, 40],
        signature_type='mismatch',
        mismatch_pattern='C->T',
        min_rate=0.05,  # Subtle signature
        rt_stop_required=False,
        confidence_weight=0.8
    ),

    'i6A': ModificationProfile(
        name='i6A',
        full_name='N6-isopentenyladenosine',
        typical_positions=[37],
        signature_type='mismatch',
        mismatch_pattern='A->G',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.0
    ),

    'm2G': ModificationProfile(
        name='m2G',
        full_name='N2-methylguanosine',
        typical_positions=[10, 26],
        signature_type='mismatch',
        mismatch_pattern='G->A',
        min_rate=0.08,
        rt_stop_required=False,
        confidence_weight=0.9
    ),

    'm22G': ModificationProfile(
        name='m22G',
        full_name='N2,N2-dimethylguanosine',
        typical_positions=[26],
        signature_type='mismatch',
        mismatch_pattern='G->A',
        min_rate=0.10,
        rt_stop_required=False,
        confidence_weight=1.0
    ),
}


class ModificationCaller:
    """
    Call specific tRNA modifications based on RT signature patterns.

    This class takes RT signature data and matches it against known modification
    profiles to identify likely modification sites with confidence scores.

    Example:
        >>> from trnaseq.modifications import RTSignatureAnalyzer, ModificationCaller
        >>>
        >>> # Get RT signatures
        >>> analyzer = RTSignatureAnalyzer()
        >>> results = analyzer.analyze_all_trnas()
        >>>
        >>> # Call modifications
        >>> caller = ModificationCaller(organism='human')
        >>> modifications = caller.call_modifications_for_trna(
        ...     'tRNA-Thr-AGT-1-1',
        ...     results['tRNA-Thr-AGT-1-1']['signatures']
        ... )
    """

    def __init__(
        self,
        organism: str = 'human',
        min_confidence: float = 0.5,
        use_position_priors: bool = True,
        statistical_test: bool = True,
        alpha: float = 0.01
    ):
        """
        Initialize modification caller.

        Args:
            organism: Organism name (for position priors)
            min_confidence: Minimum confidence score to report (0-1)
            use_position_priors: Use known modification positions to boost confidence
            statistical_test: Perform binomial test for significance
            alpha: Significance level for statistical test
        """
        self.organism = organism
        self.min_confidence = min_confidence
        self.use_position_priors = use_position_priors
        self.statistical_test = statistical_test
        self.alpha = alpha

        # Load modification profiles
        self.profiles = MODIFICATION_PROFILES

    def match_mismatch_pattern(
        self,
        pattern: str,
        ref_nt: str,
        pscm_row: pd.Series
    ) -> Tuple[bool, float]:
        """
        Check if observed mismatch pattern matches expected pattern.

        Args:
            pattern: Expected pattern (e.g., 'A->G', 'A->any', 'C->T')
            ref_nt: Reference nucleotide
            pscm_row: Row from PSCM with nucleotide counts

        Returns:
            Tuple of (matches, fraction_of_pattern)
        """
        if pattern is None:
            return True, 0.0

        parts = pattern.split('->')
        if len(parts) != 2:
            return False, 0.0

        expected_ref, expected_obs = parts

        # Check reference nucleotide matches
        if expected_ref != ref_nt:
            return False, 0.0

        # Calculate total mismatches
        total_coverage = pscm_row.sum()
        if total_coverage == 0:
            return False, 0.0

        correct_count = pscm_row.get(ref_nt, 0)
        total_mismatches = total_coverage - correct_count

        if total_mismatches == 0:
            return False, 0.0

        # Check if specific nucleotide pattern matches
        if expected_obs == 'any':
            # Any mismatch is acceptable
            fraction = total_mismatches / total_coverage
            return True, fraction
        else:
            # Specific nucleotide required
            obs_count = pscm_row.get(expected_obs, 0)
            if obs_count == 0:
                return False, 0.0
            fraction = obs_count / total_coverage
            # Pattern matches if this specific mismatch is dominant
            matches = obs_count >= (total_mismatches * 0.5)  # At least 50% of mismatches
            return matches, fraction

    def calculate_confidence(
        self,
        profile: ModificationProfile,
        position: int,
        mismatch_rate: float,
        rt_stop_pct: float,
        gap_rate: float,
        pattern_fraction: float,
        coverage: float
    ) -> float:
        """
        Calculate confidence score for modification call.

        Confidence is based on:
        1. Mismatch rate (higher is better)
        2. RT stop percentage (if required)
        3. Pattern specificity
        4. Position match (if using priors)
        5. Coverage (higher is more reliable)

        Args:
            profile: Modification profile
            position: Position in tRNA
            mismatch_rate: Observed mismatch rate
            rt_stop_pct: Observed RT stop percentage
            gap_rate: Observed gap rate
            pattern_fraction: Fraction matching expected pattern
            coverage: Read coverage

        Returns:
            Confidence score (0-1)
        """
        confidence = 0.0

        # Base confidence from mismatch rate
        if mismatch_rate >= profile.min_rate:
            # Scale between min_rate and 0.5 (saturate at 50% mismatch)
            confidence += min(1.0, (mismatch_rate - profile.min_rate) / (0.5 - profile.min_rate)) * 0.4

        # RT stop contribution
        if profile.rt_stop_required:
            if rt_stop_pct >= profile.min_rt_stop_pct:
                rt_stop_score = min(1.0, rt_stop_pct / 50.0)  # Saturate at 50%
                confidence += rt_stop_score * 0.3
            else:
                # Penalize if RT stop required but not present
                confidence *= 0.3

        # Pattern specificity
        if pattern_fraction > 0:
            confidence += pattern_fraction * 0.2

        # Position prior
        if self.use_position_priors and position in profile.typical_positions:
            confidence += 0.2

        # Coverage contribution (higher coverage = more reliable)
        coverage_score = min(1.0, np.log10(coverage + 1) / 4.0)  # Saturate at 10,000 reads
        confidence *= (0.5 + 0.5 * coverage_score)  # Scale between 50-100% based on coverage

        # Apply profile-specific weight
        confidence *= profile.confidence_weight

        # Cap at 1.0
        return min(1.0, confidence)

    def perform_statistical_test(
        self,
        coverage: int,
        mismatch_count: int,
        expected_error_rate: float = 0.01
    ) -> float:
        """
        Perform binomial test to check if mismatch rate is significantly elevated.

        Args:
            coverage: Total read coverage
            mismatch_count: Number of mismatches observed
            expected_error_rate: Expected sequencing error rate (default 1%)

        Returns:
            P-value from binomial test
        """
        if coverage == 0:
            return 1.0

        # Binomial test: is mismatch rate significantly > error rate?
        result = binomtest(
            k=int(mismatch_count),
            n=int(coverage),
            p=expected_error_rate,
            alternative='greater'
        )

        return result.pvalue

    def call_modification_at_position(
        self,
        position: int,
        signature_row: pd.Series,
        pscm_row: Optional[pd.Series] = None,
        ref_nt: Optional[str] = None
    ) -> List[Dict]:
        """
        Call modifications at a specific position.

        Args:
            position: Position in tRNA (1-based)
            signature_row: Row from signature DataFrame
            pscm_row: Row from PSCM (for pattern matching)
            ref_nt: Reference nucleotide at this position

        Returns:
            List of modification calls with confidence scores
        """
        calls = []

        mismatch_rate = signature_row.get('mismatch_rate', 0)
        rt_stop_pct = signature_row.get('rt_stop_pct', 0)
        gap_rate = signature_row.get('gap_rate', 0)
        coverage = signature_row.get('coverage', 0)

        # Try each modification profile
        for mod_name, profile in self.profiles.items():
            # Check if mismatch rate meets minimum
            if mismatch_rate < profile.min_rate:
                continue

            # Check if RT stop is required
            if profile.rt_stop_required and rt_stop_pct < profile.min_rt_stop_pct:
                continue

            # Check mismatch pattern if PSCM provided
            pattern_matches = True
            pattern_fraction = 0.0
            if pscm_row is not None and ref_nt is not None and profile.mismatch_pattern:
                pattern_matches, pattern_fraction = self.match_mismatch_pattern(
                    profile.mismatch_pattern,
                    ref_nt,
                    pscm_row
                )
                if not pattern_matches:
                    continue

            # Calculate confidence
            confidence = self.calculate_confidence(
                profile,
                position,
                mismatch_rate,
                rt_stop_pct,
                gap_rate,
                pattern_fraction,
                coverage
            )

            # Statistical test if requested
            pvalue = None
            if self.statistical_test and coverage > 0:
                mismatch_count = mismatch_rate * coverage
                pvalue = self.perform_statistical_test(coverage, mismatch_count)

                # Reduce confidence if not significant
                if pvalue > self.alpha:
                    confidence *= 0.5

            # Only report if meets minimum confidence
            if confidence >= self.min_confidence:
                call = {
                    'position': position,
                    'modification': profile.name,
                    'full_name': profile.full_name,
                    'confidence': confidence,
                    'mismatch_rate': mismatch_rate,
                    'rt_stop_pct': rt_stop_pct,
                    'gap_rate': gap_rate,
                    'coverage': coverage,
                    'pattern_fraction': pattern_fraction,
                    'pvalue': pvalue,
                    'in_typical_position': position in profile.typical_positions
                }
                calls.append(call)

        return calls

    def call_modifications_for_trna(
        self,
        trna_name: str,
        signatures_df: pd.DataFrame,
        pscm_df: Optional[pd.DataFrame] = None,
        reference_seq: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Call modifications for a single tRNA.

        Args:
            trna_name: Name of the tRNA
            signatures_df: RT signature DataFrame from RTSignatureAnalyzer
            pscm_df: Position-Specific Count Matrix (optional, for pattern matching)
            reference_seq: Reference sequence (optional, for pattern matching)

        Returns:
            DataFrame with modification calls
        """
        all_calls = []

        # Process each position with a signature
        for _, row in signatures_df.iterrows():
            if not row.get('has_signature', False):
                continue

            position = int(row['position'])

            # Get PSCM row and reference nt if available
            pscm_row = None
            ref_nt = None
            if pscm_df is not None:
                pscm_row = pscm_df.iloc[position - 1]  # Convert to 0-based
            if reference_seq is not None and position <= len(reference_seq):
                ref_nt = reference_seq[position - 1]  # Convert to 0-based

            # Call modifications at this position
            position_calls = self.call_modification_at_position(
                position, row, pscm_row, ref_nt
            )

            for call in position_calls:
                call['trna_name'] = trna_name
                all_calls.append(call)

        if not all_calls:
            return pd.DataFrame()

        # Convert to DataFrame and sort by confidence
        calls_df = pd.DataFrame(all_calls)
        calls_df = calls_df.sort_values(['confidence', 'mismatch_rate'], ascending=False)

        return calls_df

    def call_modifications_for_all_trnas(
        self,
        rt_signature_results: Dict[str, Dict[str, pd.DataFrame]],
        pscm_dict: Optional[Dict[str, pd.DataFrame]] = None,
        reference_sequences: Optional[Dict[str, str]] = None
    ) -> pd.DataFrame:
        """
        Call modifications for all tRNAs.

        Args:
            rt_signature_results: Results from RTSignatureAnalyzer.analyze_all_trnas()
            pscm_dict: Dictionary of tRNA -> PSCM DataFrame
            reference_sequences: Dictionary of tRNA -> sequence

        Returns:
            Combined DataFrame with all modification calls
        """
        all_calls = []

        for trna_name, results in rt_signature_results.items():
            signatures_df = results['signatures']

            pscm_df = pscm_dict.get(trna_name) if pscm_dict else None
            ref_seq = reference_sequences.get(trna_name) if reference_sequences else None

            trna_calls = self.call_modifications_for_trna(
                trna_name,
                signatures_df,
                pscm_df,
                ref_seq
            )

            if not trna_calls.empty:
                all_calls.append(trna_calls)

        if not all_calls:
            return pd.DataFrame()

        combined_df = pd.concat(all_calls, ignore_index=True)
        return combined_df

    def filter_by_confidence(
        self,
        calls_df: pd.DataFrame,
        min_confidence: float
    ) -> pd.DataFrame:
        """
        Filter modification calls by confidence threshold.

        Args:
            calls_df: DataFrame with modification calls
            min_confidence: Minimum confidence threshold

        Returns:
            Filtered DataFrame
        """
        return calls_df[calls_df['confidence'] >= min_confidence].copy()

    def summarize_modifications(
        self,
        calls_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Summarize modification calls by type.

        Args:
            calls_df: DataFrame with modification calls

        Returns:
            Summary DataFrame with counts per modification type
        """
        if calls_df.empty:
            return pd.DataFrame()

        summary = calls_df.groupby('modification').agg({
            'position': 'count',
            'confidence': ['mean', 'std', 'min', 'max'],
            'mismatch_rate': 'mean',
            'coverage': 'mean'
        }).round(3)

        summary.columns = ['_'.join(col).strip('_') for col in summary.columns]
        summary = summary.rename(columns={'position_count': 'num_sites'})
        summary = summary.reset_index()

        return summary
