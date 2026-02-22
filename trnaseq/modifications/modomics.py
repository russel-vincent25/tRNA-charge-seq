"""
MODOMICS Database Integration for tRNA Modification Annotation

This module integrates the MODOMICS database (https://genesilico.pl/modomics/) with
RT signature analysis to annotate observed signatures with known modification identities.

MODOMICS provides curated tRNA modification data including:
- Modified nucleotide positions in Sprinzl numbering
- Modification chemical identities and symbols
- Organism-specific modification profiles

Key workflow:
1. Fetch or load known modifications from MODOMICS for a given organism
2. Map MODOMICS Sprinzl positions to linear reference FASTA positions via alignment
3. Left-merge known modifications onto observed RT signature positions
4. Report which signatures correspond to known modifications

References:
- Boccaletto et al. 2022 (NAR) - MODOMICS database
- Sprinzl et al. 1998 (NAR) - tRNA position numbering
"""

import json
import warnings
import urllib.request
import urllib.error
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from Bio import Align


# ---------------------------------------------------------------------------
# MODOMICS symbol lookup tables
# ---------------------------------------------------------------------------

# Map MODOMICS modification symbols to their parent (unmodified) base.
# This is used to strip modification annotations from MODOMICS sequence strings
# so that the underlying base sequence can be aligned to a reference FASTA.
MODOMICS_TO_BASE: Dict[str, str] = {
    'D': 'U',       # Dihydrouridine -> U
    'T': 'U',       # m5U (ribothymidine) -> U
    'P': 'U',       # Pseudouridine (Psi) -> U
    '7': 'G',       # m7G -> G
    'K': 'G',       # m1G -> G
    '4': 'U',       # s4U -> U
    'R': 'G',       # m2,2G -> G
    'L': 'G',       # m2G -> G
    '\u0462': 'A',  # Cyrillic Yat = m1A -> A
    '\u0429': 'C',  # Cyrillic Shcha = m3C -> C
    '?': 'C',       # m5C -> C
    'I': 'A',       # Inosine -> A (deaminated)
    'Q': 'G',       # Queuosine -> G
}

# Map modification symbols to their short names.
MODOMICS_SYMBOL_NAMES: Dict[str, str] = {
    'D': 'D',        # Dihydrouridine
    'T': 'm5U',      # 5-methyluridine
    'P': 'Psi',      # Pseudouridine
    '7': 'm7G',      # 7-methylguanosine
    'K': 'm1G',      # 1-methylguanosine
    '4': 's4U',      # 4-thiouridine
    'R': 'm22G',     # N2,N2-dimethylguanosine
    'L': 'm2G',      # N2-methylguanosine
    '\u0462': 'm1A', # 1-methyladenosine
    '\u0429': 'm3C', # 3-methylcytidine
    '?': 'm5C',      # 5-methylcytidine
    'I': 'I',        # Inosine
    'Q': 'Q',        # Queuosine
}

# Map modification short names to their expected RT signature type.
# This allows prediction of what kind of RT perturbation a known modification
# should produce in sequencing data.
MOD_RT_SIGNATURES: Dict[str, str] = {
    'D': 'rt_stop',       # Dihydrouridine causes RT stops
    'm5U': 'silent',      # Ribothymidine usually silent
    'Psi': 'silent',      # Pseudouridine silent without CMC treatment
    'm7G': 'mismatch',    # G->A mismatch
    'm1G': 'rt_stop',     # Strong RT stop
    's4U': 'mismatch',    # Misincorporation
    'm22G': 'mismatch',   # G->A mismatch
    'm2G': 'mismatch',    # Subtle G->A
    'm1A': 'combined',    # Strong RT stop + A->any mismatch
    'm3C': 'mismatch',    # C->T mismatch
    'm5C': 'silent',      # Very subtle (needs bisulfite)
    'I': 'mismatch',      # A->G mismatch
    'Q': 'mismatch',      # Mismatch pattern
    'i6A': 'mismatch',    # A->G mismatch
    'ms2i6A': 'combined', # RT stop + mismatch
    't6A': 'mismatch',    # Mismatch
}

# Canonical RNA bases used to identify modification symbols.
_CANONICAL_BASES = set('ACGUacgu')

# Organism name to fallback CSV filename mapping.
_ORGANISM_CSV_MAP: Dict[str, str] = {
    'Escherichia coli': 'ecoli_known_modifications.csv',
    'Homo sapiens': 'human_known_modifications.csv',
    'Mus musculus': 'mouse_known_modifications.csv',
}

# MODOMICS API base URL.
_MODOMICS_API_BASE = 'https://genesilico.pl/modomics/api'


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _cache_dir() -> Path:
    """Return the trnaseq cache directory, creating it if necessary."""
    cache = Path.home() / '.cache' / 'trnaseq'
    cache.mkdir(parents=True, exist_ok=True)
    return cache


def _cache_path(organism: str) -> Path:
    """Return the cache file path for a given organism's MODOMICS data."""
    safe_name = organism.replace(' ', '_').lower()
    return _cache_dir() / f'modomics_{safe_name}_trna.json'


def _strip_modomics_sequence(modomics_seq: str) -> Tuple[str, List[Tuple[int, str, str]]]:
    """Strip modification symbols from a MODOMICS sequence string.

    MODOMICS encodes modified nucleotides using special Unicode characters.
    This function replaces each modification symbol with its parent base and
    records the positions and identities of the modifications.

    Args:
        modomics_seq: MODOMICS sequence string containing modification symbols.

    Returns:
        Tuple of:
        - base_seq: The plain RNA sequence (canonical bases only, U not T).
        - mods: List of (position_0based, mod_symbol, mod_short_name) tuples.
    """
    base_chars: List[str] = []
    mods: List[Tuple[int, str, str]] = []
    pos = 0

    for char in modomics_seq:
        if char in _CANONICAL_BASES:
            base_chars.append(char.upper())
            pos += 1
        elif char in MODOMICS_TO_BASE:
            parent_base = MODOMICS_TO_BASE[char]
            base_chars.append(parent_base)
            mod_name = MODOMICS_SYMBOL_NAMES.get(char, char)
            mods.append((pos, char, mod_name))
            pos += 1
        # Skip characters that are neither canonical bases nor known mod symbols
        # (e.g. whitespace or annotation characters).

    base_seq = ''.join(base_chars)
    return base_seq, mods


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class MODOMICSAnnotator:
    """Integrate MODOMICS modification database with RT signature analysis.

    This class provides methods to:
    - Fetch known tRNA modifications from the MODOMICS API (with local caching)
    - Fall back to shipped CSV files when the API is unavailable
    - Map MODOMICS Sprinzl positions to linear reference FASTA positions
    - Annotate observed RT signatures with their known modification identity

    Attributes:
        organism: Organism name (e.g. 'Escherichia coli').
        fallback_dir: Directory containing fallback CSV files.

    Example:
        >>> annotator = MODOMICSAnnotator('Escherichia coli')
        >>> mods_df = annotator.get_modifications()
        >>> annotated = annotator.annotate_signatures(signature_df, 'ecoli_tRNA-Ala-GGC-1-1')
    """

    def __init__(self, organism: str, fallback_dir: Optional[Union[str, Path]] = None):
        """
        Initialize the MODOMICS annotator.

        Args:
            organism: Organism name as it appears in MODOMICS
                      (e.g. 'Escherichia coli', 'Homo sapiens', 'Mus musculus').
            fallback_dir: Directory containing fallback CSV files. Defaults to
                          the ``data/`` subdirectory next to this module.
        """
        self.organism = organism
        self.fallback_dir = Path(fallback_dir) if fallback_dir else Path(__file__).parent / 'data'
        self._modifications_df: Optional[pd.DataFrame] = None
        self._mod_symbol_map: Optional[Dict[str, str]] = None

    # ------------------------------------------------------------------
    # Data retrieval
    # ------------------------------------------------------------------

    def fetch_modifications(self) -> pd.DataFrame:
        """Fetch tRNA modifications from the MODOMICS API.

        The method follows this resolution order:
        1. Check local disk cache (``~/.cache/trnaseq/modomics_{organism}_trna.json``).
        2. Query the MODOMICS REST API for tRNA sequences.
        3. Query the modifications endpoint for the symbol-to-name mapping.

        Successful API responses are cached locally so that subsequent calls
        do not require network access.

        Returns:
            DataFrame with columns:
                tRNA_isotype, tRNA_anticodon, position, modification_short_name,
                modification_full_name, modomics_abbrev, rt_signature_type

        Raises:
            RuntimeError: If neither cache nor API is available.
        """
        # 1. Try local cache first
        cache_file = _cache_path(self.organism)
        if cache_file.exists():
            try:
                with open(cache_file, 'r', encoding='utf-8') as fh:
                    cached = json.load(fh)
                df = pd.DataFrame(cached)
                if not df.empty:
                    self._modifications_df = df
                    return df
            except (json.JSONDecodeError, KeyError, ValueError):
                # Corrupt cache -- fall through to API
                pass

        # 2. Fetch from MODOMICS API
        organism_encoded = urllib.request.quote(self.organism)

        # Fetch tRNA sequences for this organism
        seq_url = (
            f'{_MODOMICS_API_BASE}/sequences'
            f'?RNAtype=tRNA&organism={organism_encoded}&format=json'
        )
        try:
            req = urllib.request.Request(seq_url, headers={'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                seq_data = json.loads(resp.read().decode('utf-8'))
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
            raise RuntimeError(
                f"MODOMICS API unavailable for organism '{self.organism}': {exc}"
            ) from exc

        # 3. Fetch symbol -> full name mapping
        mod_map = self._fetch_modification_names()

        # Parse sequences into modification records
        records: List[Dict] = []
        entries = seq_data if isinstance(seq_data, list) else seq_data.get('results', [])

        for entry in entries:
            modomics_seq = entry.get('sequence', entry.get('seq', ''))
            subtype = entry.get('subtype', '')
            anticodon = entry.get('anticodon', '')

            # Parse isotype from subtype (e.g. "tRNA-Ala" -> "Ala")
            isotype = subtype.replace('tRNA-', '').replace('tRNA', '').strip()
            if not isotype:
                isotype = entry.get('isotype', 'unknown')

            # Strip modification symbols to find positions
            _base_seq, mods = _strip_modomics_sequence(modomics_seq)

            for pos_0, _symbol, short_name in mods:
                full_name = mod_map.get(short_name, short_name)
                rt_sig = MOD_RT_SIGNATURES.get(short_name, 'unknown')

                records.append({
                    'tRNA_isotype': isotype,
                    'tRNA_anticodon': anticodon,
                    'position': pos_0 + 1,  # 1-based
                    'modification_short_name': short_name,
                    'modification_full_name': full_name,
                    'modomics_abbrev': _symbol,
                    'rt_signature_type': rt_sig,
                })

        df = pd.DataFrame(records)

        # Cache results
        try:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_file, 'w', encoding='utf-8') as fh:
                json.dump(records, fh, ensure_ascii=False, indent=2)
        except OSError:
            pass  # Non-fatal -- caching is best-effort

        self._modifications_df = df
        return df

    def _fetch_modification_names(self) -> Dict[str, str]:
        """Fetch modification short-name to full-name mapping from MODOMICS API.

        Returns:
            Dictionary mapping short names (e.g. 'm1A') to full names
            (e.g. '1-methyladenosine'). Falls back to a built-in mapping
            if the API call fails.
        """
        if self._mod_symbol_map is not None:
            return self._mod_symbol_map

        # Built-in fallback mapping
        builtin: Dict[str, str] = {
            'D': 'dihydrouridine',
            'm5U': '5-methyluridine (ribothymidine)',
            'Psi': 'pseudouridine',
            'm7G': '7-methylguanosine',
            'm1G': '1-methylguanosine',
            's4U': '4-thiouridine',
            'm22G': 'N2,N2-dimethylguanosine',
            'm2G': 'N2-methylguanosine',
            'm1A': '1-methyladenosine',
            'm3C': '3-methylcytidine',
            'm5C': '5-methylcytidine',
            'I': 'inosine',
            'Q': 'queuosine',
            'i6A': 'N6-isopentenyladenosine',
            'ms2i6A': '2-methylthio-N6-isopentenyladenosine',
            't6A': 'N6-threonylcarbamoyladenosine',
        }

        mod_url = f'{_MODOMICS_API_BASE}/modifications?format=json'
        try:
            req = urllib.request.Request(mod_url, headers={'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                mod_data = json.loads(resp.read().decode('utf-8'))

            entries = mod_data if isinstance(mod_data, list) else mod_data.get('results', [])
            for entry in entries:
                short = entry.get('short_name', entry.get('abbreviation', ''))
                full = entry.get('name', entry.get('full_name', ''))
                if short and full:
                    builtin[short] = full
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError):
            # API unavailable -- use built-in mapping only
            pass

        self._mod_symbol_map = builtin
        return builtin

    def load_fallback(self, organism: Optional[str] = None) -> pd.DataFrame:
        """Load modification data from a shipped CSV fallback file.

        This is used when the MODOMICS API is unavailable. The CSV files are
        stored in the ``data/`` subdirectory next to this module.

        Args:
            organism: Organism name override. Uses ``self.organism`` if None.

        Returns:
            DataFrame with the same schema as :meth:`fetch_modifications`.

        Raises:
            FileNotFoundError: If no fallback CSV exists for the organism.
        """
        org = organism or self.organism

        csv_name = _ORGANISM_CSV_MAP.get(org)
        if csv_name is None:
            # Try a sanitised version of the organism name
            safe = org.replace(' ', '_').lower()
            csv_name = f'{safe}_known_modifications.csv'

        csv_path = self.fallback_dir / csv_name
        if not csv_path.exists():
            raise FileNotFoundError(
                f"No fallback CSV found for organism '{org}' at {csv_path}. "
                f"Available organisms: {list(_ORGANISM_CSV_MAP.keys())}"
            )

        df = pd.read_csv(csv_path, keep_default_na=False)

        # Normalise column names to the canonical schema
        rename_map = {}
        expected_cols = {
            'tRNA_isotype', 'tRNA_anticodon', 'position',
            'modification_short_name', 'modification_full_name',
            'modomics_abbrev', 'rt_signature_type',
        }
        for col in df.columns:
            lower = col.lower().strip()
            for expected in expected_cols:
                if lower == expected.lower():
                    rename_map[col] = expected
        if rename_map:
            df = df.rename(columns=rename_map)

        # Ensure rt_signature_type is populated from MOD_RT_SIGNATURES
        if 'rt_signature_type' not in df.columns:
            df['rt_signature_type'] = df['modification_short_name'].map(
                MOD_RT_SIGNATURES
            ).fillna('unknown')

        self._modifications_df = df
        return df

    def get_modifications(self, use_api: bool = True) -> pd.DataFrame:
        """Get tRNA modifications, trying the API first then the fallback.

        This is the recommended entry point for obtaining modification data.
        It provides a seamless experience: fast when cached, resilient when
        the network is unavailable.

        Args:
            use_api: If True, try the MODOMICS API (with caching) before
                     falling back to CSV. If False, use the fallback directly.

        Returns:
            DataFrame with columns: tRNA_isotype, tRNA_anticodon, position,
            modification_short_name, modification_full_name, modomics_abbrev,
            rt_signature_type.
        """
        if self._modifications_df is not None:
            return self._modifications_df

        if use_api:
            try:
                return self.fetch_modifications()
            except (RuntimeError, Exception) as exc:
                warnings.warn(
                    f"MODOMICS API fetch failed ({exc}); falling back to CSV.",
                    stacklevel=2,
                )

        # Fallback to shipped CSV
        try:
            return self.load_fallback()
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"No modification data available for '{self.organism}'. "
                f"API failed and no fallback CSV exists. Original error: {exc}"
            ) from exc

    # ------------------------------------------------------------------
    # Position mapping
    # ------------------------------------------------------------------

    def map_to_reference(
        self,
        modomics_seq: str,
        ref_seq: str
    ) -> List[Tuple[int, int, str, str]]:
        """Map MODOMICS sequence positions to linear reference FASTA positions.

        MODOMICS sequences use Sprinzl numbering, which does not correspond
        directly to linear FASTA positions due to variable-length loops.
        This method aligns the MODOMICS base sequence to the reference FASTA
        and maps modification positions through the alignment coordinates.

        Algorithm:
        1. Strip modification symbols from the MODOMICS sequence to recover
           the underlying base sequence (replacing each mod symbol with its
           parent base via :data:`MODOMICS_TO_BASE`).
        2. Align the base sequence to the reference using
           ``Bio.Align.PairwiseAligner`` (global mode, DNA-friendly scoring).
        3. Walk the best alignment to translate each MODOMICS position that
           carries a modification into the corresponding linear reference
           position.

        Args:
            modomics_seq: MODOMICS sequence string with modification symbols.
            ref_seq: Reference FASTA sequence (DNA or RNA, will be uppercased).

        Returns:
            List of tuples ``(modomics_pos_1based, linear_pos_1based,
            mod_short_name, rt_signature_type)`` for every modification found
            in the MODOMICS sequence. Positions that fall in alignment gaps
            (i.e. insertions in the MODOMICS sequence relative to the reference)
            are assigned a linear position of -1.
        """
        # Step 1: strip modification symbols
        base_seq, mods = _strip_modomics_sequence(modomics_seq)

        if not mods:
            return []

        # Convert ref_seq to RNA uppercase for consistent comparison
        ref_upper = ref_seq.upper().replace('T', 'U')
        base_seq_upper = base_seq.upper()

        # Step 2: align
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1

        alignments = aligner.align(ref_upper, base_seq_upper)
        if len(alignments) == 0:
            warnings.warn(
                "Could not align MODOMICS sequence to reference; "
                "returning unmapped positions.",
                stacklevel=2,
            )
            results = []
            for pos_0, _symbol, short_name in mods:
                rt_sig = MOD_RT_SIGNATURES.get(short_name, 'unknown')
                results.append((pos_0 + 1, -1, short_name, rt_sig))
            return results

        best = alignments[0]

        # Step 3: build query-pos -> target-pos mapping from alignment blocks
        target_blocks, query_blocks = best.aligned
        query_to_target: Dict[int, int] = {}
        for (t_start, t_end), (q_start, q_end) in zip(target_blocks, query_blocks):
            for offset in range(t_end - t_start):
                query_to_target[q_start + offset] = t_start + offset

        # Map modification positions
        results: List[Tuple[int, int, str, str]] = []
        for pos_0, _symbol, short_name in mods:
            rt_sig = MOD_RT_SIGNATURES.get(short_name, 'unknown')
            linear_pos_0 = query_to_target.get(pos_0, -2)
            # Convert to 1-based; -1 sentinel for unmapped
            linear_pos_1 = linear_pos_0 + 1 if linear_pos_0 >= 0 else -1
            results.append((pos_0 + 1, linear_pos_1, short_name, rt_sig))

        return results

    # ------------------------------------------------------------------
    # Annotation
    # ------------------------------------------------------------------

    def annotate_signatures(
        self,
        signature_df: pd.DataFrame,
        trna_name: str
    ) -> pd.DataFrame:
        """Merge observed RT signatures with known MODOMICS modifications.

        Performs a left join so that every position in the input signature
        DataFrame is preserved. Positions that coincide with known
        modifications receive annotation columns; others get NaN.

        Args:
            signature_df: DataFrame from
                :meth:`~trnaseq.modifications.rt_signatures.RTSignatureAnalyzer.identify_signature_positions`.
                Must contain a ``position`` column (1-based linear positions).
            trna_name: tRNA name used to look up known modifications. The name
                is parsed using the convention
                ``{prefix}_tRNA-{aa}-{anticodon}-{copy}-{allele}``.

        Returns:
            A copy of *signature_df* with added columns:

            - ``known_modification``: Short name of the known modification at
              this position (or NaN).
            - ``modification_full_name``: Full chemical name.
            - ``expected_rt_signature``: Expected RT signature type from
              MODOMICS (e.g. 'mismatch', 'rt_stop', 'combined', 'silent').
            - ``is_known_modification``: Boolean flag.
        """
        known_df = self.get_known_modifications(trna_name)

        result = signature_df.copy()

        if known_df.empty:
            result['known_modification'] = np.nan
            result['modification_full_name'] = np.nan
            result['expected_rt_signature'] = np.nan
            result['is_known_modification'] = False
            return result

        # Prepare known modifications for merge (keep relevant columns)
        merge_cols = ['position', 'modification_short_name',
                      'modification_full_name', 'rt_signature_type']
        known_subset = known_df[[c for c in merge_cols if c in known_df.columns]].copy()

        # Deduplicate: if multiple mods at same position, concatenate names
        if known_subset.duplicated(subset='position').any():
            grouped = known_subset.groupby('position').agg({
                'modification_short_name': lambda x: '; '.join(sorted(set(x))),
                'modification_full_name': lambda x: '; '.join(sorted(set(x)))
                    if 'modification_full_name' in known_subset.columns else '',
                'rt_signature_type': lambda x: '; '.join(sorted(set(x)))
                    if 'rt_signature_type' in known_subset.columns else '',
            }).reset_index()
            known_subset = grouped

        # Rename for merge
        rename = {
            'modification_short_name': 'known_modification',
            'modification_full_name': 'modification_full_name',
            'rt_signature_type': 'expected_rt_signature',
        }
        known_subset = known_subset.rename(columns=rename)

        # Ensure position dtype matches
        result['position'] = result['position'].astype(int)
        known_subset['position'] = known_subset['position'].astype(int)

        # Left merge
        result = result.merge(known_subset, on='position', how='left')

        # Add boolean flag
        result['is_known_modification'] = result['known_modification'].notna()

        return result

    def get_known_modifications(self, trna_name: str) -> pd.DataFrame:
        """Get all known modifications for a specific tRNA.

        Parses the tRNA name to extract the amino acid isotype and anticodon,
        then filters the loaded modification database accordingly.

        Name format convention::

            {prefix}_tRNA-{aa}-{anticodon}-{copy}-{allele}

        Parsing: ``aa = name.split('-')[1]``, ``anticodon = name.split('-')[2]``.

        Args:
            trna_name: tRNA name following the project naming convention.

        Returns:
            DataFrame of known modifications for this tRNA with their
            reference positions. Returns an empty DataFrame if the tRNA
            isotype/anticodon combination is not found in the database.
        """
        # Parse tRNA name
        parts = trna_name.split('-')
        if len(parts) < 3:
            warnings.warn(
                f"Cannot parse tRNA name '{trna_name}' — expected at least "
                f"3 hyphen-separated parts (got {len(parts)}). "
                f"Returning empty modifications.",
                stacklevel=2,
            )
            return pd.DataFrame()

        isotype = parts[1]     # e.g. 'Ala'
        anticodon = parts[2]   # e.g. 'GGC'

        # Ensure modifications are loaded
        mods_df = self.get_modifications()
        if mods_df.empty:
            return pd.DataFrame()

        # Filter by isotype and anticodon (case-insensitive matching)
        mask = (
            mods_df['tRNA_isotype'].str.lower() == isotype.lower()
        ) & (
            mods_df['tRNA_anticodon'].str.upper() == anticodon.upper()
        )
        matched = mods_df[mask].copy()

        if matched.empty:
            # Try matching by isotype only (some databases lack anticodon)
            mask_iso = mods_df['tRNA_isotype'].str.lower() == isotype.lower()
            matched = mods_df[mask_iso].copy()

        return matched.reset_index(drop=True)

    # ------------------------------------------------------------------
    # Representations
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        n_mods = len(self._modifications_df) if self._modifications_df is not None else 0
        status = f'{n_mods} modifications loaded' if n_mods else 'not loaded'
        return f"MODOMICSAnnotator(organism='{self.organism}', {status})"
