# RVHMS04_5215_library

The **5,215 tRNA library as designed and ordered** — i.e. the sequences physically
present in the plasmid pool that was cloned and sequenced.

|  |  |
|---|---|
| sequences | 5,216 = 5,164 tdb + 3 synthetic spike-in controls + 49 *E. coli* K12 natives |
| length | 67–96 nt |
| md5 (`.fa`) | `5341afc9c67656b374e57206cb2be16d` |
| built | 2026-08-31 with `makeblastdb -dbtype prot -in RVHMS04_5215_library.fa -blastdb_version 4` |
| BLASTDB version | 4 (required by SWIPE) |

## Why this exists — do not use `RVHMS04_5215` for this data

`RVHMS04_5215.fa` and this file have **identical headers, identical order, and the
same 52 controls**. They differ in **276 of 5,164** tdb sequences.

`RVHMS04_5215.fa` carries in-silico corrections that were **never synthesised**:

* a G−1 added to 120 His tRNAs (correct for *E. coli* HisRS, but absent from the DNA —
  *E. coli* has no Thg1, so the transcripts cannot acquire it)
* genomic 3′ ends preserved instead of the design's standardised CCA

This file instead carries the sequences that were actually ordered, verified three ways:

1. byte-identical to the 2023-11-16 amplicon reference (5,164/5,164)
2. byte-identical to the surviving TWIST oPool order files
   (2023-12-11 subsets 4,553/4,553 · 2024-01-18 CCA 2,110/2,110 · 2024-08-30 116/116)
3. confirmed in the raw reads — the `+CCCA` tail present, the His G−1 absent

Aligning against `RVHMS04_5215` mismatches the molecules at those 276 positions.

## Provenance

Generated in the RVHMS04-NovelOTS repo:

    python scripts/update_library_fasta.py --library 5215 --out data/processed/refs

Per-sequence before/after: `data/processed/refs/RVHMS04_5215_library_changes.tsv`.
The matching proK-flanked amplicon reference (for amplicon-seq, not tRNA-seq) is
`RVHMS04_5215_amplicons.fa` in the same repo.

`RVHMS04_5215_redundant.csv` is carried over unchanged from `RVHMS04_5215/` —
membership is identical, only sequences changed.
