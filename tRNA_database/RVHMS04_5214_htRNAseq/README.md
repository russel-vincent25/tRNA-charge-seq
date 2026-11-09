# RVHMS04_5214_htRNAseq

`RVHMS04_5215_library` with the two synthetic controls that were never spiked into the
2023-11-07 h-tRNAseq removed. Use this for that experiment.

|  |  |
|---|---|
| sequences | 5,214 = 5,164 tdbD + 49 *E. coli* K12 natives + 1 synthetic spike-in |
| length | 67–96 nt |
| md5 (`.fa`) | `17655cc6b5c5cafa89df861faa2944ec` |
| built | 2026-11-09 with `makeblastdb -dbtype prot -in RVHMS04_5214_htRNAseq.fa -blastdb_version 4` |
| BLASTDB version | 4 (required by SWIPE) |

## What changed from RVHMS04_5215_library

Exactly two records removed:

    Synthetic_Control_tRNA-Phe-AAA-1-1
    Synthetic_Control_tRNA-Val-TAC-1-1

Every retained sequence is byte-identical, and headers and record order are unchanged.

`Synthetic_Control_tRNA-Tyr-CTA-1-1` is **retained and load-bearing** — it is the
Mj-Tyr-CUA deacylation spike-in that serves as the per-sample assay-floor control.
RusselLib lacked it, which is why its reads were absorbed by native `tdbD00003510`
(3,779/3,425/4,083 UMI at ~0% charge). Do not remove it.

## Why the two were dropped

Phe-AAA and Val-TAC were controls for the ±aaRS experiment, not this one. In
`results/charge/charge_df_transcript.csv` from the 5215_library run they have **zero
rows in both arms**, while Tyr-CTA carries 79,319 UMI (ssR1.33) and 107,950 (ssR1.35).
They are not in this library, so they only offered two extra landing sites for a
mismapped read.

## Inherited verification

Because the retained sequences are untouched, this file inherits
`RVHMS04_5215_library`'s verification. It carries the **as-ordered** sequences, not the
in-silico corrections in `RVHMS04_5215.fa` (which added a G−1 to His tRNAs and kept
genomic 3′ ends instead of the standardised CCA). That was confirmed three ways: against
the 2023-11-16 amplicon reference (5,164/5,164), against the surviving TWIST oPool order
files, and in the raw reads.

**Do not use `RVHMS04_5215` or `RusselLib` for this data.** Both disagree with the
synthesised DNA in 276 of 5,164 sequences, 146 of them at the 3′ terminus, which is
exactly where charge is read.

## Provenance

Built in the RVHMS04-NovelOTS repo from
`internal/data/processed/refs/RVHMS04_5215_library.fa`; the copy here is
`internal/data/processed/refs/RVHMS04_5214_htRNAseq.fa`.

Matching pipeline config and sample list:

    internal/exp/2023-11-07-h-tRNASeq/config_h-tRNAseq_5214_htRNAseq.yaml
    internal/exp/2023-11-07-h-tRNASeq/data/sample_list_RVHMS04B-htRNASeq_5214_htRNAseq.xlsx

`RVHMS04_5214_htRNAseq_redundant.csv` is carried over unchanged from
`RVHMS04_5215_library/` — construct membership is identical.
