"""
Single-read Analysis of Crosstalks (SLAC) for tRNA Modifications

Detects coordination between modifications at distinct positions within
the same tRNA molecule by analyzing per-read mismatch patterns.

For each pair of known modification positions on a tRNA, builds a 2x2
contingency table from individual reads:

                    Position B
                  Modified  Unmodified
    Position A  ┌─────────┬───────────┐
    Modified    │   n11   │    n10    │
    Unmodified  │   n01   │    n00    │
                └─────────┴───────────┘

Fisher's exact test yields an odds ratio and p-value:
- OR > 1: modifications tend to co-occur (positive coordination)
- OR < 1: modifications are anti-correlated (negative coordination)
- OR ≈ 1: modifications are independent

References:
- Behrens et al. 2023 (NAR) — SLAC methodology in mim-tRNAseq
"""

import bz2
import json
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict
from scipy.stats import fisher_exact


def _is_mismatch(qseq_char: str, dseq_char: str) -> bool:
    """Check if a position shows a mismatch (potential modification signal)."""
    if qseq_char == '-' or dseq_char == '-':
        return False  # Indels handled separately
    return qseq_char.upper() != dseq_char.upper()


def _extract_read_mismatches(
    qseq: str,
    dseq: str,
    dpos_start: int,
    positions_of_interest: Set[int],
) -> Dict[int, bool]:
    """Extract mismatch status at specific positions from a single read.

    Args:
        qseq: Aligned query (read) sequence.
        dseq: Aligned reference sequence.
        dpos_start: 1-based start position on reference.
        positions_of_interest: Set of 1-based reference positions to check.

    Returns:
        Dict mapping each covered position of interest to True (mismatch)
        or False (match). Positions not covered by the read are omitted.
    """
    result = {}
    ref_pos = dpos_start  # 1-based

    for i in range(len(qseq)):
        if dseq[i] == '-':
            # Insertion in read — ref position doesn't advance
            continue
        if ref_pos in positions_of_interest:
            result[ref_pos] = _is_mismatch(qseq[i], dseq[i])
        ref_pos += 1

    return result


class CrosstalkAnalyzer:
    """Analyze modification crosstalks from per-read alignment data.

    This class streams SWalign JSON files and builds contingency tables
    for all pairs of known modification positions on each tRNA.

    Attributes:
        min_coverage: Minimum reads covering both positions in a pair.
        mismatch_threshold: Minimum mismatch rate to consider a position
            as potentially modified (filters out noise positions).
        alpha: Significance threshold for Fisher's exact test.
    """

    def __init__(
        self,
        min_coverage: int = 50,
        mismatch_threshold: float = 0.05,
        alpha: float = 0.05,
    ):
        self.min_coverage = min_coverage
        self.mismatch_threshold = mismatch_threshold
        self.alpha = alpha

    def analyze_sample(
        self,
        json_path: Path,
        mod_positions: Dict[str, List[int]],
        min_fmax_score: float = 0.8,
    ) -> pd.DataFrame:
        """Analyze modification crosstalks for one sample.

        Args:
            json_path: Path to SWalign JSON (optionally bz2-compressed).
            mod_positions: Dict mapping tRNA name to list of 1-based
                modification positions to test. Typically from
                ``MODOMICSAnnotator.get_known_mods_linear()``.
            min_fmax_score: Minimum fractional alignment score to include
                a read.

        Returns:
            DataFrame with columns: trna_name, pos_a, pos_b, mod_a, mod_b,
            n_both, n_a_only, n_b_only, n_neither, odds_ratio, pvalue,
            log2_odds_ratio, significant, coordination_type.
        """
        # Build contingency tables per tRNA per position pair
        # contingency[(trna, pos_a, pos_b)] = [n_both, n_a_only, n_b_only, n_neither]
        contingency: Dict[Tuple[str, int, int], List[int]] = defaultdict(
            lambda: [0, 0, 0, 0]
        )

        # Stream the JSON
        open_fn = bz2.open if str(json_path).endswith('.bz2') else open
        try:
            import json_stream
            use_stream = True
        except ImportError:
            use_stream = False

        if use_stream:
            with open_fn(json_path, 'rt', encoding='utf-8') as fh:
                data = json_stream.load(fh)
                for read_id, align_dict in data.persistent().items():
                    self._process_read(
                        align_dict, mod_positions, min_fmax_score, contingency
                    )
        else:
            with open_fn(json_path, 'rt', encoding='utf-8') as fh:
                data = json.load(fh)
            for read_id, align_dict in data.items():
                self._process_read(
                    align_dict, mod_positions, min_fmax_score, contingency
                )

        # Build results DataFrame
        return self._build_results(contingency, mod_positions)

    def _process_read(
        self,
        align_dict: dict,
        mod_positions: Dict[str, List[int]],
        min_fmax_score: float,
        contingency: Dict[Tuple[str, int, int], List[int]],
    ) -> None:
        """Process a single read's alignment data."""
        # Filter: must be aligned with sufficient quality
        if not align_dict.get('aligned', False):
            return
        if align_dict.get('Fmax_score', 0) < min_fmax_score:
            return

        # Get tRNA name (first match if multi-mapped)
        trna_name = str(align_dict.get('name', ''))
        if '@' in trna_name:
            trna_name = trna_name.split('@')[0]

        # Check if we have modification positions for this tRNA
        positions = mod_positions.get(trna_name)
        if not positions or len(positions) < 2:
            return

        qseq = str(align_dict.get('qseq', ''))
        dseq = str(align_dict.get('dseq', ''))
        dpos_raw = align_dict.get('dpos', [0, 0])
        # json_stream returns PersistentStreamingJSONList; convert to plain list
        dpos = list(dpos_raw) if not isinstance(dpos_raw, (list, tuple)) else dpos_raw
        dpos_start = int(dpos[0])

        if not qseq or not dseq or dpos_start < 1:
            return

        positions_set = set(positions)

        # Get mismatch status at each position of interest
        mm_status = _extract_read_mismatches(qseq, dseq, dpos_start, positions_set)

        # For each pair of positions both covered by this read, update contingency
        covered = sorted(mm_status.keys())
        for i in range(len(covered)):
            for j in range(i + 1, len(covered)):
                pos_a, pos_b = covered[i], covered[j]
                a_mod = mm_status[pos_a]
                b_mod = mm_status[pos_b]

                key = (trna_name, pos_a, pos_b)
                if a_mod and b_mod:
                    contingency[key][0] += 1  # both modified
                elif a_mod and not b_mod:
                    contingency[key][1] += 1  # A only
                elif not a_mod and b_mod:
                    contingency[key][2] += 1  # B only
                else:
                    contingency[key][3] += 1  # neither

    def _build_results(
        self,
        contingency: Dict[Tuple[str, int, int], List[int]],
        mod_positions: Dict[str, List[int]],
    ) -> pd.DataFrame:
        """Convert contingency tables to a results DataFrame."""
        rows = []

        for (trna_name, pos_a, pos_b), counts in contingency.items():
            n_both, n_a_only, n_b_only, n_neither = counts
            total = sum(counts)

            if total < self.min_coverage:
                continue

            # Filter: both positions must show some modification signal
            rate_a = (n_both + n_a_only) / total if total > 0 else 0
            rate_b = (n_both + n_b_only) / total if total > 0 else 0
            if rate_a < self.mismatch_threshold or rate_b < self.mismatch_threshold:
                continue

            # Fisher's exact test
            table = np.array([[n_both, n_a_only], [n_b_only, n_neither]])
            try:
                odds_ratio, pvalue = fisher_exact(table)
            except ValueError:
                odds_ratio, pvalue = np.nan, np.nan

            # Log2 odds ratio (with pseudocount to avoid log(0))
            if odds_ratio > 0 and np.isfinite(odds_ratio):
                log2_or = np.log2(odds_ratio)
            else:
                log2_or = np.nan

            # Classify coordination type
            if pvalue > self.alpha or np.isnan(pvalue):
                coord_type = 'independent'
            elif odds_ratio > 1:
                coord_type = 'positive'  # co-occurring
            else:
                coord_type = 'negative'  # anti-correlated

            rows.append({
                'trna_name': trna_name,
                'pos_a': pos_a,
                'pos_b': pos_b,
                'rate_a': round(rate_a, 4),
                'rate_b': round(rate_b, 4),
                'n_both': n_both,
                'n_a_only': n_a_only,
                'n_b_only': n_b_only,
                'n_neither': n_neither,
                'n_total': total,
                'odds_ratio': round(odds_ratio, 4) if np.isfinite(odds_ratio) else np.nan,
                'log2_odds_ratio': round(log2_or, 4) if np.isfinite(log2_or) else np.nan,
                'pvalue': pvalue,
                'significant': pvalue < self.alpha if np.isfinite(pvalue) else False,
                'coordination_type': coord_type,
            })

        if not rows:
            return pd.DataFrame()

        df = pd.DataFrame(rows)

        # FDR correction (Benjamini-Hochberg)
        if len(df) > 1:
            pvals = df['pvalue'].values
            valid = np.isfinite(pvals)
            if valid.any():
                from trnaseq.modifications.modification_caller import benjamini_hochberg_fdr
                fdr_sig = np.full(len(pvals), False)
                fdr_sig[valid] = benjamini_hochberg_fdr(pvals[valid], alpha=self.alpha)
                df['fdr_significant'] = fdr_sig
            else:
                df['fdr_significant'] = False
        else:
            df['fdr_significant'] = df['significant']

        return df.sort_values('pvalue').reset_index(drop=True)

    def analyze_multiple_samples(
        self,
        json_paths: Dict[str, Path],
        mod_positions: Dict[str, List[int]],
        min_fmax_score: float = 0.8,
    ) -> pd.DataFrame:
        """Analyze crosstalks across multiple samples and aggregate.

        Args:
            json_paths: Dict mapping sample name to SWalign JSON path.
            mod_positions: Dict mapping tRNA name to modification positions.
            min_fmax_score: Minimum fractional alignment score.

        Returns:
            DataFrame with per-sample and aggregated crosstalk results.
        """
        all_results = []

        for sample_name, json_path in json_paths.items():
            try:
                result = self.analyze_sample(
                    json_path, mod_positions, min_fmax_score
                )
                if not result.empty:
                    result['sample'] = sample_name
                    all_results.append(result)
            except Exception as exc:
                warnings.warn(
                    f"Crosstalk analysis failed for {sample_name}: {exc}",
                    stacklevel=2,
                )

        if not all_results:
            return pd.DataFrame()

        return pd.concat(all_results, ignore_index=True)
