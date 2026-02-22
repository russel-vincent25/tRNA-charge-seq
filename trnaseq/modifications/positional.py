"""
Position-Specific Count Matrix (PSCM) Extraction

Streams SWalign JSON files and builds per-position nucleotide count matrices
for each tRNA reference. Supports RT stop profiling, mismatch analysis,
and anticodon coverage assessment.

PSCM array layout per tRNA -- ndarray of shape (ref_len, 8):
    Col 0: A count
    Col 1: C count
    Col 2: G count
    Col 3: T count
    Col 4: N count
    Col 5: gap (deletion) count
    Col 6: coverage (total non-insertion observations)
    Col 7: rt_stop count (reads whose 5' end maps here)
"""

import numpy as np
import pandas as pd
import bz2
import json_stream
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Union
from Bio import SeqIO
from mpire import WorkerPool

# Column indices in the PSCM array
_COL = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4,
        'gap': 5, 'coverage': 6, 'rt_stop': 7}
_NT_COLS = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
_NCOLS = 8

# All 12 possible single-nucleotide substitutions
_SUBSTITUTIONS = [
    f'{r}_to_{a}' for r in 'ACGT' for a in 'ACGT' if r != a
]


def _load_reference(fasta_path: Union[str, Path]) -> Dict[str, dict]:
    """Load reference tRNA sequences from a FASTA file.

    Returns:
        {name: {'seq': str, 'seq_len': int}}
    """
    ref_dict = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        ref_dict[record.id] = {
            'seq': str(record.seq),
            'seq_len': len(record.seq),
        }
    return ref_dict


def _parse_anticodon_from_name(trna_name: str) -> Optional[str]:
    """Extract anticodon from tRNA name using the standard naming convention.

    Name format: {prefix}_tRNA-{aa}-{anticodon}-{copy}-{allele}
    Parser uses split('-')[2] for the anticodon.
    """
    parts = trna_name.split('-')
    if len(parts) >= 3:
        return parts[2]
    return None


class PositionalExtractor:
    """Extract per-position count matrices from SWalign JSON files.

    Each aligned read is walked character-by-character against the reference
    to build a Position-Specific Count Matrix (PSCM) tracking nucleotide
    observations, deletions, coverage, and RT stop positions.

    Example::

        ext = PositionalExtractor('ecoli_tRNAs.fa', min_fmax_score=0.8)
        pscm = ext.extract_sample('sample1_SWalign.json.bz2')
        mismatch_df = ext.compute_mismatch_profile(pscm)
        rt_df = ext.compute_rt_profile(pscm)
    """

    def __init__(
        self,
        reference_fasta: Union[str, Path],
        min_fmax_score: float = 0.0,
        only_unique: bool = True,
    ):
        """
        Args:
            reference_fasta: Path to tRNA FASTA file.
            min_fmax_score: Minimum Fmax_score filter (0.0 = no filter).
            only_unique: If True, skip multi-mapped reads (name contains '@').
        """
        self.reference_fasta = Path(reference_fasta)
        self.min_fmax_score = min_fmax_score
        self.only_unique = only_unique
        self.ref_dict = _load_reference(self.reference_fasta)

    # ------------------------------------------------------------------
    # Core extraction
    # ------------------------------------------------------------------
    def extract_sample(self, json_path: Union[str, Path]) -> Dict[str, np.ndarray]:
        """Stream a SWalign JSON and build a PSCM per tRNA.

        Args:
            json_path: Path to ``{sample}_SWalign.json.bz2``.

        Returns:
            ``{tRNA_name: ndarray(ref_len, 8)}`` -- see module docstring
            for column layout.

        Algorithm per aligned read:
            1. Filter on aligned flag, uniqueness, and Fmax_score.
            2. Walk ``qseq`` / ``dseq`` char-by-char starting at ``dpos[0]``.
            3. Insertion (``dseq[i]=='-'``): skip, ref_pos stays.
            4. Deletion (``qseq[i]=='-'``): increment gap + coverage, advance ref_pos.
            5. Match/mismatch: increment nucleotide + coverage, advance ref_pos.
            6. Record RT stop at ``dpos[0]`` (5' end of read on reference).
        """
        json_path = Path(json_path)
        pscm: Dict[str, np.ndarray] = {}

        with bz2.open(json_path, 'rt', encoding='utf-8') as fh:
            data = json_stream.load(fh)
            for _read_id, align_dict in data.persistent().items():
                if not align_dict['aligned']:
                    continue

                name: str = align_dict['name']

                # Multi-mapped filter
                if self.only_unique and '@' in name:
                    continue

                # Fmax score filter
                fmax = float(align_dict['Fmax_score'])
                if fmax < self.min_fmax_score:
                    continue

                # Resolve tRNA name (first before '@' if multi-mapped)
                trna_name = name.split('@')[0]
                if trna_name not in self.ref_dict:
                    warnings.warn(
                        f"tRNA '{trna_name}' from alignment not found in reference; skipping read."
                    )
                    continue

                ref_len = self.ref_dict[trna_name]['seq_len']

                # Lazily initialise PSCM
                if trna_name not in pscm:
                    pscm[trna_name] = np.zeros((ref_len, _NCOLS), dtype=np.float64)
                mat = pscm[trna_name]

                qseq: str = align_dict['qseq']
                dseq: str = align_dict['dseq']
                dpos = align_dict['dpos']  # [start, end], 1-indexed
                ref_pos = int(dpos[0])  # 1-indexed

                # RT stop: 5' boundary on reference
                idx_rt = ref_pos - 1
                if 0 <= idx_rt < ref_len:
                    mat[idx_rt, _COL['rt_stop']] += 1

                # Walk alignment character-by-character
                for i in range(len(qseq)):
                    d_char = dseq[i]
                    q_char = qseq[i]

                    if d_char == '-':
                        # Insertion in read vs reference -- skip
                        continue

                    # Map 1-indexed ref_pos to 0-indexed array row
                    idx = ref_pos - 1
                    if 0 <= idx < ref_len:
                        if q_char == '-':
                            # Deletion in read
                            mat[idx, _COL['gap']] += 1
                        else:
                            # Match or mismatch
                            nt_col = _NT_COLS.get(q_char.upper())
                            if nt_col is not None:
                                mat[idx, nt_col] += 1
                            else:
                                # Treat any unexpected character as N
                                mat[idx, _COL['N']] += 1
                        # Coverage for all non-insertion cases
                        mat[idx, _COL['coverage']] += 1

                    ref_pos += 1

        return pscm

    # ------------------------------------------------------------------
    # RT-stop profile
    # ------------------------------------------------------------------
    def compute_rt_profile(
        self,
        pscm_dict: Dict[str, np.ndarray],
    ) -> pd.DataFrame:
        """Per-tRNA, per-position RT stop profile.

        Returns a DataFrame with columns:
            tRNA_name, position (1-based), rt_stop_count, rt_stop_fraction,
            cumulative_coverage, ref_nt

        ``rt_stop_fraction`` is computed as rt_stop_count divided by the
        coverage at that position. ``cumulative_coverage`` is the running
        sum of coverage from the 3' end to each position (i.e. total reads
        that passed through or stopped at that position).
        """
        rows: List[dict] = []
        for trna_name, mat in pscm_dict.items():
            ref_seq = self.ref_dict[trna_name]['seq']
            ref_len = mat.shape[0]

            coverage = mat[:, _COL['coverage']]
            rt_stop = mat[:, _COL['rt_stop']]

            # Cumulative coverage from 3' end (reverse cumsum)
            cum_cov = np.cumsum(coverage[::-1])[::-1]

            for pos_idx in range(ref_len):
                cov = coverage[pos_idx]
                rt_cnt = rt_stop[pos_idx]
                rt_frac = rt_cnt / cov if cov > 0 else 0.0
                rows.append({
                    'tRNA_name': trna_name,
                    'position': pos_idx + 1,
                    'rt_stop_count': int(rt_cnt),
                    'rt_stop_fraction': rt_frac,
                    'cumulative_coverage': int(cum_cov[pos_idx]),
                    'ref_nt': ref_seq[pos_idx] if pos_idx < len(ref_seq) else 'N',
                })

        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Mismatch profile
    # ------------------------------------------------------------------
    def compute_mismatch_profile(
        self,
        pscm_dict: Dict[str, np.ndarray],
    ) -> pd.DataFrame:
        """Per-tRNA, per-position mismatch and deletion rates.

        Returns a DataFrame with columns:
            tRNA_name, position (1-based), ref_nt, coverage,
            mismatch_rate, deletion_rate,
            A_count, C_count, G_count, T_count, N_count, gap_count,
            plus 12 substitution columns (A_to_C, A_to_G, ..., T_to_G).

        ``mismatch_rate = (coverage - correct_nt_count - gap_count) / coverage``

        Positions where the reference nucleotide is 'N' are treated as
        all-match (mismatch_rate = 0).
        """
        rows: List[dict] = []
        for trna_name, mat in pscm_dict.items():
            ref_seq = self.ref_dict[trna_name]['seq']
            ref_len = mat.shape[0]

            for pos_idx in range(ref_len):
                cov = mat[pos_idx, _COL['coverage']]
                a_cnt = mat[pos_idx, _COL['A']]
                c_cnt = mat[pos_idx, _COL['C']]
                g_cnt = mat[pos_idx, _COL['G']]
                t_cnt = mat[pos_idx, _COL['T']]
                n_cnt = mat[pos_idx, _COL['N']]
                gap_cnt = mat[pos_idx, _COL['gap']]

                ref_nt = ref_seq[pos_idx] if pos_idx < len(ref_seq) else 'N'

                # Mismatch rate
                if cov > 0 and ref_nt.upper() in _NT_COLS:
                    correct_cnt = mat[pos_idx, _NT_COLS[ref_nt.upper()]]
                    mm_rate = (cov - correct_cnt - gap_cnt) / cov
                else:
                    mm_rate = 0.0

                # N-masked reference: treat as all-match
                if ref_nt.upper() == 'N':
                    mm_rate = 0.0

                del_rate = gap_cnt / cov if cov > 0 else 0.0

                row = {
                    'tRNA_name': trna_name,
                    'position': pos_idx + 1,
                    'ref_nt': ref_nt,
                    'coverage': int(cov),
                    'mismatch_rate': mm_rate,
                    'deletion_rate': del_rate,
                    'A_count': int(a_cnt),
                    'C_count': int(c_cnt),
                    'G_count': int(g_cnt),
                    'T_count': int(t_cnt),
                    'N_count': int(n_cnt),
                    'gap_count': int(gap_cnt),
                }

                # Per-substitution counts (e.g. A_to_C)
                nt_counts = {'A': a_cnt, 'C': c_cnt, 'G': g_cnt, 'T': t_cnt}
                for sub in _SUBSTITUTIONS:
                    ref_b, _, alt_b = sub.split('_')
                    if ref_nt.upper() == ref_b:
                        row[sub] = int(nt_counts[alt_b])
                    else:
                        row[sub] = 0

                rows.append(row)

        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Anticodon coverage
    # ------------------------------------------------------------------
    def compute_anticodon_coverage(
        self,
        pscm_dict: Dict[str, np.ndarray],
        anticodon_positions: Optional[Dict[str, tuple]] = None,
    ) -> pd.DataFrame:
        """Fraction of reads that cover the anticodon positions.

        Args:
            anticodon_positions: ``{tRNA_name: (pos1, pos2, pos3)}`` (1-based).
                If None, auto-detect by searching for the anticodon subsequence
                in the reference sequence.

        Returns a DataFrame with columns:
            tRNA_name, anticodon_pos_start, anticodon_pos_end,
            reads_covering_anticodon, total_reads, fraction_covering,
            mean_anticodon_coverage
        """
        if anticodon_positions is None:
            anticodon_positions = self._autodetect_anticodon_positions(
                list(pscm_dict.keys())
            )

        rows: List[dict] = []
        for trna_name, mat in pscm_dict.items():
            ac_pos = anticodon_positions.get(trna_name)
            if ac_pos is None:
                warnings.warn(
                    f"No anticodon positions for '{trna_name}'; skipping."
                )
                continue

            # Convert 1-based positions to 0-based indices
            indices = [p - 1 for p in ac_pos]
            valid = [i for i in indices if 0 <= i < mat.shape[0]]
            if not valid:
                continue

            ac_coverages = np.array([mat[i, _COL['coverage']] for i in valid])
            mean_ac_cov = float(ac_coverages.mean())
            # Reads covering anticodon = minimum coverage across the 3 positions
            reads_covering = int(ac_coverages.min())

            # Total reads for this tRNA: max coverage at any position (proxy)
            total_reads = int(mat[:, _COL['coverage']].max())
            frac = reads_covering / total_reads if total_reads > 0 else 0.0

            rows.append({
                'tRNA_name': trna_name,
                'anticodon_pos_start': min(ac_pos),
                'anticodon_pos_end': max(ac_pos),
                'reads_covering_anticodon': reads_covering,
                'total_reads': total_reads,
                'fraction_covering': frac,
                'mean_anticodon_coverage': mean_ac_cov,
            })

        return pd.DataFrame(rows)

    def _autodetect_anticodon_positions(
        self,
        trna_names: List[str],
    ) -> Dict[str, tuple]:
        """Auto-detect anticodon positions by finding the anticodon subsequence
        in each reference sequence.

        The anticodon is extracted from the tRNA name (``split('-')[2]``).
        We search for that trinucleotide in the reference and take the first
        occurrence in the plausible range (positions 25-45 in linear coords).
        """
        result: Dict[str, tuple] = {}
        for name in trna_names:
            anticodon = _parse_anticodon_from_name(name)
            if anticodon is None or name not in self.ref_dict:
                continue

            ref_seq = self.ref_dict[name]['seq'].upper()
            ac_upper = anticodon.upper()

            # Search within a plausible range (positions 25-45, 0-indexed 24-44)
            search_start = max(0, 24)
            search_end = min(len(ref_seq), 45)
            search_region = ref_seq[search_start:search_end]
            idx = search_region.find(ac_upper)

            if idx >= 0:
                # Convert back to 1-based full-sequence positions
                pos_start = search_start + idx + 1
                result[name] = (pos_start, pos_start + 1, pos_start + 2)
            else:
                # Fallback: search entire sequence
                idx_full = ref_seq.find(ac_upper)
                if idx_full >= 0:
                    pos_start = idx_full + 1
                    result[name] = (pos_start, pos_start + 1, pos_start + 2)
                else:
                    warnings.warn(
                        f"Anticodon '{anticodon}' not found in reference for '{name}'."
                    )

        return result

    # ------------------------------------------------------------------
    # Parallel execution
    # ------------------------------------------------------------------
    def run_parallel(
        self,
        json_dir: Union[str, Path],
        sample_names: List[str],
        n_jobs: int = 4,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        """Run extraction over multiple samples in parallel using mpire.

        Args:
            json_dir: Directory containing ``{sample}_SWalign.json.bz2`` files.
            sample_names: List of sample names.
            n_jobs: Number of parallel workers.

        Returns:
            ``{sample_name: {tRNA_name: ndarray}}``
        """
        json_dir = Path(json_dir)
        json_paths = []
        valid_names = []
        for name in sample_names:
            p = json_dir / f'{name}_SWalign.json.bz2'
            if p.exists():
                json_paths.append(p)
                valid_names.append(name)
            else:
                warnings.warn(f"JSON file not found for sample '{name}': {p}")

        if not json_paths:
            warnings.warn("No valid JSON files found; returning empty dict.")
            return {}

        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self.extract_sample, json_paths)

        return {name: res for name, res in zip(valid_names, results)}
