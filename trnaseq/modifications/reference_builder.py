"""
Build tRNA Reference FASTAs from MODOMICS Isodecoder Sequences

Instead of using genomic tRNA databases (GtRNAdb/tRNAscan-SE) which contain
redundant gene copies, this module builds reference FASTAs from MODOMICS
isodecoder sequences — unique mature tRNA transcripts with embedded
modification annotations.

Benefits:
- Eliminates multi-mapping from gene copy redundancy
- Each sequence is a unique isodecoder (same anticodon, distinct body)
- Modification positions are inherently known from the MODOMICS symbols
- Works for any organism in MODOMICS (E. coli, human, mouse, yeast, etc.)

Usage:
    >>> from trnaseq.modifications.reference_builder import build_modomics_reference
    >>> ref_path = build_modomics_reference('human', output_dir='/path/to/refs')
    # Produces: /path/to/refs/human_modomics.fa
    #           /path/to/refs/human_modomics_modifications.csv
"""

import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from .modomics import (
    MODOMICSAnnotator,
    _strip_modomics_sequence,
    MODOMICS_SYMBOL_NAMES,
    MODOMICS_TO_BASE,
)


def build_modomics_reference(
    organism: str,
    output_dir: Union[str, Path],
    use_api: bool = True,
    name_prefix: Optional[str] = None,
) -> Tuple[Path, Path]:
    """Build a reference FASTA and modification table from MODOMICS sequences.

    For each unique isodecoder in MODOMICS, produces:
    - A DNA FASTA entry (RNA U→T converted)
    - A row in the modification CSV for each modified position

    Args:
        organism: Organism name (e.g. 'ecoli', 'human', 'mouse').
        output_dir: Directory to write output files.
        use_api: Whether to try the MODOMICS API (True) or use CSV only.
        name_prefix: Optional prefix for FASTA headers. Defaults to organism.

    Returns:
        Tuple of (fasta_path, mod_csv_path).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    annotator = MODOMICSAnnotator(organism)
    annotator.get_modifications(use_api=use_api)

    if not annotator._modomics_sequences:
        raise ValueError(
            f"No MODOMICS sequences available for '{organism}'. "
            f"Ensure the MODOMICS API is accessible or provide sequences manually."
        )

    prefix = name_prefix or annotator.organism.replace(' ', '_')
    safe_org = organism.lower().replace(' ', '_')

    fasta_path = output_dir / f'{safe_org}_modomics.fa'
    mod_csv_path = output_dir / f'{safe_org}_modomics_modifications.csv'

    fasta_entries: List[str] = []
    mod_rows: List[Dict] = []

    for (isotype, anticodon), modomics_seq in sorted(annotator._modomics_sequences.items()):
        # Strip modifications to get base sequence
        base_seq_rna, mods = _strip_modomics_sequence(modomics_seq)
        base_seq_dna = base_seq_rna.replace('U', 'T').replace('u', 't')

        # Build FASTA header: prefix_tRNA-Iso-ACN-1-1
        iso_cap = isotype.capitalize()
        trna_name = f'{prefix}_tRNA-{iso_cap}-{anticodon}-1-1'

        fasta_entries.append(f'>{trna_name}\n{base_seq_dna}')

        # Record modifications
        for pos_0, symbol, short_name in mods:
            parent_base = MODOMICS_TO_BASE.get(symbol, 'N')
            identity = base_seq_dna[pos_0] if pos_0 < len(base_seq_dna) else 'N'

            mod_rows.append({
                'tRNA': trna_name,
                'pos': pos_0 + 1,  # 1-based
                'identity': identity,
                'base': parent_base,
                'mod': symbol,
                'mod_name': short_name,
            })

    # Write FASTA
    with open(fasta_path, 'w') as fh:
        fh.write('\n'.join(fasta_entries) + '\n')

    # Write modification CSV (same format as all_tRNA_pos_mod_info.csv)
    with open(mod_csv_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['tRNA', 'pos', 'identity', 'base', 'mod', 'mod_name'])
        writer.writeheader()
        writer.writerows(mod_rows)

    return fasta_path, mod_csv_path


def build_full_position_table(
    organism: str,
    output_path: Union[str, Path],
    use_api: bool = True,
    name_prefix: Optional[str] = None,
) -> Path:
    """Build a complete per-position modification table (like all_tRNA_pos_mod_info.csv).

    For every position in every MODOMICS isodecoder, records:
    - tRNA name, position, FASTA identity (DNA), RNA base, modification symbol

    Unmodified positions show the canonical base in the 'mod' column.

    Args:
        organism: Organism name.
        output_path: Path for the output CSV.
        use_api: Whether to try the MODOMICS API.
        name_prefix: Optional prefix for tRNA names.

    Returns:
        Path to the written CSV.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    annotator = MODOMICSAnnotator(organism)
    annotator.get_modifications(use_api=use_api)

    if not annotator._modomics_sequences:
        raise ValueError(f"No MODOMICS sequences for '{organism}'.")

    prefix = name_prefix or annotator.organism.replace(' ', '_')
    rows: List[Dict] = []

    for (isotype, anticodon), modomics_seq in sorted(annotator._modomics_sequences.items()):
        base_seq_rna, mods = _strip_modomics_sequence(modomics_seq)
        base_seq_dna = base_seq_rna.replace('U', 'T').replace('u', 't')

        iso_cap = isotype.capitalize()
        trna_name = f'{prefix}_tRNA-{iso_cap}-{anticodon}-1-1'

        # Build a mod lookup: pos_0 → (symbol, short_name)
        mod_at = {}
        for pos_0, symbol, short_name in mods:
            mod_at[pos_0] = symbol

        for pos_0 in range(len(base_seq_dna)):
            identity = base_seq_dna[pos_0]
            base_rna = base_seq_rna[pos_0]
            mod_symbol = mod_at.get(pos_0, base_rna)

            rows.append({
                'tRNA': trna_name,
                'pos': pos_0 + 1,
                'identity': identity,
                'base': base_rna,
                'mod': mod_symbol,
            })

    with open(output_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['tRNA', 'pos', 'identity', 'base', 'mod'])
        writer.writeheader()
        writer.writerows(rows)

    return output_path
