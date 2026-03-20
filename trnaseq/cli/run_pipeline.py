#!/usr/bin/env python3
"""
Unified Preprocessing Pipeline (Stages 0-7)
============================================

Combines all preprocessing and analysis steps into a single script:
- Stage 0a: Adapter Removal & Merging
- Stage 0b: Barcode Splitting
- Stage 0c: UMI Trimming
- Stage 1: SWIPE Alignment
- Stage 2: Stats Collection
- Stage 3: Charge Quantification (optional) + Fragment Analysis
- Stage 4: Parquet Storage (optional)
- Stage 5: QC Summary Report (optional)
- Stage 6: Modification Analysis (optional)
- Stage 7: Differential Abundance Analysis (optional)

Output:
- inp_file_df.xlsx (input file summary)
- sample_df.xlsx (sample information with QC metrics)
- logs/pipeline.log, logs/computing_metrics.csv (operational files)
- data/stats_collection/ALL_stats_aggregate.csv (combined statistics)
- results/charge/ (charge quantification + report, if enabled)
- results/fragments/ (fragment classification + RT drop-off + report, if enabled)
- results/parquet/ (parquet-format data, if enabled)
- results/modifications/ (modification calls + report, if enabled)
- results/abundance/ (differential abundance + report, if enabled)
- qc_reports/ (QC dashboard, if enabled)

Usage:
    python -m trnaseq pipeline \
        --config config.yaml \
        --project-dir /path/to/project/ \
        --n-jobs 8 \
        --parquet

Author: Based on projects/example/process_data.ipynb
Date: 2026-02-11
"""

import os
import sys
import time
import json
import shutil
from enum import Enum
from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
import argparse
import yaml
import pandas as pd
import numpy as np

# Add repo root to path for src imports
repo_path = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_path))

from src.misc import (
    index_to_sample_df,
    downsample_raw_input,
    read_tRNAdb_info,
    sample_df_to_dict
)
from src.read_processing import AR_merge, BC_split, UMI_trim
from src.alignment import SWIPE_align
from src.stats_collection import STATS_collection

# Optional imports for extended analysis stages
try:
    from trnaseq.charge.quantifier import ChargeQuantifier
    CHARGE_AVAILABLE = True
except ImportError:
    CHARGE_AVAILABLE = False

try:
    from trnaseq.io.storage import tRNAseqDataStore, PYARROW_AVAILABLE
    STORAGE_AVAILABLE = PYARROW_AVAILABLE
except ImportError:
    STORAGE_AVAILABLE = False

try:
    from trnaseq.qc.report import QCReportGenerator
    QC_AVAILABLE = True
except ImportError:
    QC_AVAILABLE = False

try:
    from trnaseq.fragments import FragmentAnalyser
    FRAGMENTS_AVAILABLE = True
except ImportError:
    FRAGMENTS_AVAILABLE = False

try:
    from trnaseq.modifications.positional import PositionalExtractor
    from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer
    from trnaseq.modifications.modification_caller import (
        ModificationCaller,
        estimate_background_error_rate,
        ReplicateAggregator,
    )
    from trnaseq.modifications.modomics import MODOMICSAnnotator
    MODIFICATIONS_AVAILABLE = True
except ImportError:
    MODIFICATIONS_AVAILABLE = False

try:
    from trnaseq.qc.modification_report import ModificationReportGenerator
    MOD_REPORT_AVAILABLE = True
except ImportError:
    MOD_REPORT_AVAILABLE = False

try:
    from trnaseq.qc.charge_report import ChargeReportGenerator
    CHARGE_REPORT_AVAILABLE = True
except ImportError:
    CHARGE_REPORT_AVAILABLE = False

try:
    from trnaseq.qc.fragment_report import FragmentReportGenerator
    FRAGMENT_REPORT_AVAILABLE = True
except ImportError:
    FRAGMENT_REPORT_AVAILABLE = False

try:
    from trnaseq.abundance import DifferentialAbundance
    ABUNDANCE_AVAILABLE = True
except ImportError:
    ABUNDANCE_AVAILABLE = False

try:
    from trnaseq.qc.abundance_report import AbundanceReportGenerator
    ABUNDANCE_REPORT_AVAILABLE = True
except ImportError:
    ABUNDANCE_REPORT_AVAILABLE = False


# ------------------------------------------------------------------
# Preflight validation helpers
# ------------------------------------------------------------------

class CheckStatus(Enum):
    PASS = "PASS"
    WARN = "WARN"
    FAIL = "FAIL"
    SKIP = "SKIP"


@dataclass
class CheckResult:
    name: str
    status: CheckStatus
    message: str = ""
    group: str = ""


class PreflightReport:
    """Collects check results and prints a grouped, color-coded report."""

    _COLORS = {
        CheckStatus.PASS: "\033[32m",  # green
        CheckStatus.WARN: "\033[33m",  # yellow
        CheckStatus.FAIL: "\033[31m",  # red
        CheckStatus.SKIP: "\033[90m",  # grey
    }
    _RESET = "\033[0m"

    def __init__(self):
        self.results: list[CheckResult] = []

    def add(self, result: CheckResult):
        self.results.append(result)

    @property
    def ok(self) -> bool:
        return not any(r.status == CheckStatus.FAIL for r in self.results)

    def print(self):
        use_color = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()
        W = 70

        def _tag(status: CheckStatus) -> str:
            label = status.value.center(4)
            if use_color:
                return f"{self._COLORS[status]}{label}{self._RESET}"
            return label

        print("=" * W)
        print("  tRNA-charge-seq PREFLIGHT CHECK".center(W))
        print("=" * W)

        groups: dict[str, list[CheckResult]] = {}
        for r in self.results:
            groups.setdefault(r.group, []).append(r)

        for group, checks in groups.items():
            print(f"\n  [{group}]")
            for c in checks:
                tag = _tag(c.status)
                # Pad name to 32 chars for alignment
                name_part = c.name[:32].ljust(32)
                info = f"  ({c.message})" if c.message and c.status != CheckStatus.FAIL else ""
                print(f"    {tag}  {name_part}{info}")
                if c.status == CheckStatus.FAIL and c.message:
                    for line in c.message.split("\n"):
                        print(f"          {line}")

        n_pass = sum(1 for r in self.results if r.status == CheckStatus.PASS)
        n_warn = sum(1 for r in self.results if r.status == CheckStatus.WARN)
        n_fail = sum(1 for r in self.results if r.status == CheckStatus.FAIL)
        n_skip = sum(1 for r in self.results if r.status == CheckStatus.SKIP)

        print()
        print("-" * W)
        verdict = "READY" if self.ok else "NOT READY"
        parts = [f"{n_pass} passed"]
        if n_warn:
            parts.append(f"{n_warn} warnings")
        if n_fail:
            parts.append(f"{n_fail} failures")
        if n_skip:
            parts.append(f"{n_skip} skipped")
        summary = ", ".join(parts)
        if use_color:
            color = self._COLORS[CheckStatus.PASS] if self.ok else self._COLORS[CheckStatus.FAIL]
            print(f"  {color}{verdict}{self._RESET}  {summary}")
        else:
            print(f"  {verdict}  {summary}")
        print("=" * W)


# ------------------------------------------------------------------
# Path fields in config that should be resolved relative to project_dir
# ------------------------------------------------------------------
_PATH_FIELDS = ['sample_list', 'index_list', 'SWIPE_score_mat', 'common_seqs', 'adapter_sequences']
_DICT_PATH_FIELDS = ['tRNA_database']  # dict values are paths


def _save_df(df, path_stem, write_csv=False):
    """Save DataFrame as parquet with CSV fallback."""
    try:
        df.to_parquet(f'{path_stem}.parquet', index=False)
    except Exception:
        df.to_csv(f'{path_stem}.csv', index=False)
        return
    if write_csv:
        df.to_csv(f'{path_stem}.csv', index=False)


def _pool_pscm_dicts(all_pscm):
    """Sum PSCM arrays across all samples for background estimation.

    Args:
        all_pscm: {sample_name: {trna_name: ndarray(ref_len, 8)}}

    Returns:
        {trna_name: ndarray(ref_len, 8)} — element-wise sum across samples.
    """
    pooled = {}
    for sample_pscm in all_pscm.values():
        for trna_name, mat in sample_pscm.items():
            if trna_name not in pooled:
                pooled[trna_name] = mat.copy()
            else:
                if pooled[trna_name].shape == mat.shape:
                    pooled[trna_name] += mat
    return pooled


class PreprocessingPipeline:
    """
    Unified preprocessing pipeline for tRNA-charge-seq

    Runs stages 0a-5 and outputs key files for downstream analysis.
    """

    def __init__(self, config_file, project_dir, n_jobs=4, sample_index=None,
                 threads_per_job=None):
        """
        Initialize pipeline

        Parameters:
            config_file: YAML configuration file
            project_dir: Project directory (must contain data/{seq_dir}/ for stages 0-1)
            n_jobs: Number of parallel jobs
            sample_index: Process only sample at this 0-based index (for SLURM array jobs)
            threads_per_job: Threads per subprocess (AR/SWIPE). Default: from config or 2.
        """
        self.config_file = Path(config_file).resolve()
        self.project_dir = Path(project_dir).resolve()
        self.n_jobs = n_jobs
        self.sample_index = sample_index
        # Create output directory
        self.project_dir.mkdir(parents=True, exist_ok=True)

        # Log file goes in logs/ subdirectory
        logs_dir = self.project_dir / 'logs'
        logs_dir.mkdir(exist_ok=True)
        self.log_file = logs_dir / "pipeline.log"

        # Load configuration
        with open(self.config_file, 'r') as f:
            self.config = yaml.safe_load(f)

        # Resolve relative paths in config against project_dir
        self._resolve_config_paths()

        # Resolve threads_per_job: CLI arg > config > default 2
        if threads_per_job is not None:
            self.threads_per_job = threads_per_job
        else:
            self.threads_per_job = self.config.get('threads_per_job', 2)

        # Initialize status tracking
        self.status = {
            'start_time': datetime.now(),
            'stages_completed': []
        }

        # Computing metrics
        self.stage_timings = {}
        self.file_metrics = {}

    def _resolve_config_paths(self):
        """Resolve relative paths in config against project_dir.

        A path is treated as relative if it doesn't start with '/'.
        Absolute paths pass through unchanged. None/null values are skipped.
        """
        for key in _PATH_FIELDS:
            val = self.config.get(key)
            if val is not None and not os.path.isabs(val):
                self.config[key] = str(self.project_dir / val)

        for key in _DICT_PATH_FIELDS:
            d = self.config.get(key)
            if isinstance(d, dict):
                for k, v in d.items():
                    if v is not None and not os.path.isabs(v):
                        d[k] = str(self.project_dir / v)

    # ------------------------------------------------------------------
    # Preflight validation
    # ------------------------------------------------------------------

    def preflight(self, stages=None):
        """Run preflight checks and print a report.

        Returns:
            0 if all checks pass, 1 if any check fails.
        """
        if stages is None:
            stages = {'0a', '0b', '0c', '1', '2', '3', '4', '5', '6', '7'}

        report = PreflightReport()

        self._preflight_config(report, stages)
        self._preflight_files(report, stages)
        self._preflight_database(report, stages)
        self._preflight_samples(report, stages)
        self._preflight_tools(report, stages)

        report.print()
        return 0 if report.ok else 1

    def _preflight_config(self, report, stages):
        """Validate configuration values."""
        group = "Config"

        # Config YAML already loaded successfully by __init__
        report.add(CheckResult("Config YAML is valid", CheckStatus.PASS, group=group))

        # umi_trim_mode
        umi_mode = self.config.get('umi_trim_mode', 'anchored')
        valid_modes = ('anchored', 'pyrimidine')
        if umi_mode in valid_modes:
            report.add(CheckResult(
                "umi_trim_mode", CheckStatus.PASS,
                f"{umi_mode}", group=group))
        else:
            report.add(CheckResult(
                "umi_trim_mode", CheckStatus.FAIL,
                f"'{umi_mode}' is not valid. Use: {', '.join(valid_modes)}",
                group=group))

        # umi_anchor when anchored
        if umi_mode == 'anchored':
            anchor = self.config.get('umi_anchor')
            if anchor:
                # Check if it exists in adapter_sequences
                adapter_path = self.config.get('adapter_sequences')
                if adapter_path and Path(adapter_path).is_file():
                    with open(adapter_path) as f:
                        adapter_seqs = yaml.safe_load(f)
                    if anchor in adapter_seqs:
                        report.add(CheckResult(
                            "umi_anchor", CheckStatus.PASS,
                            f"{anchor} = {adapter_seqs[anchor]}", group=group))
                    else:
                        report.add(CheckResult(
                            "umi_anchor", CheckStatus.FAIL,
                            f"'{anchor}' not in {adapter_path}\n"
                            f"Available: {list(adapter_seqs.keys())}",
                            group=group))
                else:
                    report.add(CheckResult(
                        "umi_anchor", CheckStatus.PASS,
                        f"{anchor} (adapter_sequences file checked later)",
                        group=group))
            else:
                report.add(CheckResult(
                    "umi_anchor", CheckStatus.FAIL,
                    "'umi_anchor' required when umi_trim_mode='anchored'",
                    group=group))

        # seq_dir: warn if contains 'data/' (double-nested)
        seq_dir = self.config.get('seq_dir', 'raw_fastq')
        if 'data/' in seq_dir or 'data\\' in seq_dir:
            report.add(CheckResult(
                "seq_dir", CheckStatus.WARN,
                f"'{seq_dir}' contains 'data/' — may be double-nested.\n"
                f"Pipeline uses project_dir/data/{{seq_dir}}/",
                group=group))

        # Stage 6: needs organism
        if '6' in stages and self.config.get('run_modification_analysis', False):
            org = self.config.get('organism')
            if not org:
                report.add(CheckResult(
                    "organism (stage 6)", CheckStatus.FAIL,
                    "Stage 6 (modifications) requires 'organism' in config",
                    group=group))
            else:
                report.add(CheckResult(
                    "organism (stage 6)", CheckStatus.PASS,
                    org, group=group))

        # Stage 7: needs abundance_control
        if '7' in stages and self.config.get('run_abundance_analysis', False):
            ctrl = self.config.get('abundance_control')
            if not ctrl:
                report.add(CheckResult(
                    "abundance_control (stage 7)", CheckStatus.FAIL,
                    "Stage 7 (abundance) requires 'abundance_control' in config",
                    group=group))
            else:
                report.add(CheckResult(
                    "abundance_control (stage 7)", CheckStatus.PASS,
                    ctrl, group=group))

    def _preflight_files(self, report, stages):
        """Check that required input files exist."""
        group = "Files"

        needs_early = stages.intersection({'0a', '0b', '0c', '1'})

        # sample_list
        sample_list = self.config.get('sample_list')
        if sample_list:
            p = Path(sample_list)
            if p.is_file():
                try:
                    df = pd.read_excel(p)
                    report.add(CheckResult(
                        "sample_list", CheckStatus.PASS,
                        f"{len(df)} rows", group=group))
                except Exception as e:
                    report.add(CheckResult(
                        "sample_list", CheckStatus.FAIL,
                        f"Cannot read: {e}", group=group))
            else:
                report.add(CheckResult(
                    "sample_list", CheckStatus.FAIL,
                    f"Not found: {p}", group=group))
        else:
            report.add(CheckResult(
                "sample_list", CheckStatus.FAIL,
                "Not set in config", group=group))

        # index_list
        index_list = self.config.get('index_list')
        if index_list:
            p = Path(index_list)
            if p.is_file():
                report.add(CheckResult("index_list", CheckStatus.PASS, group=group))
            else:
                report.add(CheckResult(
                    "index_list", CheckStatus.FAIL,
                    f"Not found: {p}", group=group))

        # SWIPE_score_mat
        score_mat = self.config.get('SWIPE_score_mat')
        if score_mat and stages.intersection({'1'}):
            p = Path(score_mat)
            if p.is_file():
                try:
                    from src.alignment import read_scoremat
                    read_scoremat(str(p))
                    report.add(CheckResult(
                        "SWIPE_score_mat", CheckStatus.PASS,
                        "parsed OK", group=group))
                except Exception as e:
                    report.add(CheckResult(
                        "SWIPE_score_mat", CheckStatus.FAIL,
                        f"Parse error: {e}", group=group))
            else:
                report.add(CheckResult(
                    "SWIPE_score_mat", CheckStatus.FAIL,
                    f"Not found: {p}", group=group))

        # adapter_sequences (if anchored mode)
        umi_mode = self.config.get('umi_trim_mode', 'anchored')
        if umi_mode == 'anchored':
            adapter_path = self.config.get('adapter_sequences')
            if adapter_path:
                p = Path(adapter_path)
                if p.is_file():
                    report.add(CheckResult(
                        "adapter_sequences", CheckStatus.PASS, group=group))
                else:
                    report.add(CheckResult(
                        "adapter_sequences", CheckStatus.FAIL,
                        f"Not found: {p}", group=group))
            else:
                # Check default location
                default = self.project_dir / 'utils' / 'adapter_sequences.yaml'
                if default.is_file():
                    report.add(CheckResult(
                        "adapter_sequences", CheckStatus.PASS,
                        f"using default: {default}", group=group))
                else:
                    report.add(CheckResult(
                        "adapter_sequences", CheckStatus.FAIL,
                        f"Not set and default not found: {default}",
                        group=group))

        # common_seqs (.bz2)
        common_seqs = self.config.get('common_seqs')
        if common_seqs is not None:
            p = Path(common_seqs)
            if p.is_file():
                report.add(CheckResult("common_seqs", CheckStatus.PASS, group=group))
            else:
                report.add(CheckResult(
                    "common_seqs", CheckStatus.FAIL,
                    f"Not found: {p}", group=group))

        # Raw data directory + FASTQ files
        if needs_early:
            seq_dir = self.config.get('seq_dir', 'raw_fastq')
            raw_path = self.project_dir / 'data' / seq_dir
            if raw_path.is_dir():
                report.add(CheckResult(
                    "Raw data directory", CheckStatus.PASS,
                    str(raw_path), group=group))

                # Check FASTQ files
                try:
                    sl = self.config.get('sample_list')
                    if sl and Path(sl).is_file():
                        sdf = pd.read_excel(sl)
                        fastq_cols = [c for c in sdf.columns
                                      if 'fastq' in c.lower() and 'filename' in c.lower()]
                        missing = []
                        for col in fastq_cols:
                            for fname in sdf[col].dropna().unique():
                                fpath = raw_path / fname
                                if not fpath.is_file():
                                    missing.append(str(fname))
                                elif fpath.stat().st_size == 0:
                                    missing.append(f"{fname} (empty)")
                        if missing:
                            n = len(missing)
                            shown = missing[:5]
                            msg = f"{n} missing:\n" + "\n".join(shown)
                            if n > 5:
                                msg += f"\n... and {n - 5} more"
                            report.add(CheckResult(
                                "FASTQ files", CheckStatus.FAIL, msg, group=group))
                        else:
                            n_files = sum(
                                len(sdf[col].dropna().unique()) for col in fastq_cols)
                            report.add(CheckResult(
                                "FASTQ files", CheckStatus.PASS,
                                f"{n_files} files found", group=group))
                except Exception:
                    pass  # sample_list checks handle parse errors
            else:
                report.add(CheckResult(
                    "Raw data directory", CheckStatus.FAIL,
                    f"Not found: {raw_path}", group=group))

    def _preflight_database(self, report, stages):
        """Validate tRNA database files."""
        group = "Database"
        tRNA_db = self.config.get('tRNA_database')
        if not isinstance(tRNA_db, dict) or not tRNA_db:
            report.add(CheckResult(
                "tRNA_database", CheckStatus.FAIL,
                "Must be a dict of species: fasta_path", group=group))
            return

        for species, fasta_path in tRNA_db.items():
            p = Path(fasta_path)

            # FASTA file
            if not p.is_file():
                report.add(CheckResult(
                    f"tRNA FASTA ({species})", CheckStatus.FAIL,
                    f"Not found: {p}", group=group))
                continue

            # Count sequences and validate headers
            # Format: {prefix}_tRNA-{aa}-{anticodon}-{copy}-{allele}
            # Parser uses split('-')[1]=aa, split('-')[2]=anticodon
            # Prefix must not contain hyphens (would shift split indices)
            n_seqs = 0
            bad_headers = []
            hyphen_prefix = []
            with open(p) as fh:
                for line in fh:
                    if not line.startswith('>'):
                        continue
                    n_seqs += 1
                    header = line[1:].strip().split()[0]  # first word
                    if '_tRNA-' not in header:
                        bad_headers.append(header)
                    else:
                        prefix = header.split('_tRNA-')[0]
                        if '-' in prefix:
                            hyphen_prefix.append(header)

            report.add(CheckResult(
                f"tRNA FASTA ({species})", CheckStatus.PASS,
                f"{n_seqs} sequences", group=group))

            if bad_headers:
                shown = bad_headers[:3]
                msg = f"{len(bad_headers)} headers missing 'tRNA' segment:\n"
                msg += "\n".join(f"  {h}" for h in shown)
                report.add(CheckResult(
                    f"FASTA headers ({species})", CheckStatus.WARN,
                    msg, group=group))

            if hyphen_prefix:
                shown = hyphen_prefix[:3]
                msg = (f"{len(hyphen_prefix)} headers have hyphens in prefix "
                       f"(breaks parser):\n")
                msg += "\n".join(f"  {h}" for h in shown)
                if len(hyphen_prefix) > 3:
                    msg += f"\n  ... and {len(hyphen_prefix) - 3} more"
                report.add(CheckResult(
                    f"FASTA prefix ({species})", CheckStatus.FAIL,
                    msg, group=group))

            # BLAST protein database files (appended to full .fa path)
            blast_exts = ['.phr', '.pin', '.psq']
            missing_blast = [ext for ext in blast_exts
                             if not Path(str(p) + ext).is_file()]
            if not missing_blast:
                report.add(CheckResult(
                    f"BLAST db ({species})", CheckStatus.PASS, group=group))
            else:
                msg = f"Missing: {', '.join(missing_blast)}\n"
                msg += f"Run: makeblastdb -dbtype prot -in {p.name}"
                report.add(CheckResult(
                    f"BLAST db ({species})", CheckStatus.FAIL, msg, group=group))

    def _preflight_samples(self, report, stages):
        """Validate sample and index list contents."""
        group = "Samples"

        sample_list = self.config.get('sample_list')
        index_list = self.config.get('index_list')

        # Need both files to exist for these checks
        if not sample_list or not Path(sample_list).is_file():
            report.add(CheckResult(
                "Sample validation", CheckStatus.SKIP,
                "sample_list not available", group=group))
            return

        try:
            sdf = pd.read_excel(sample_list)
        except Exception:
            report.add(CheckResult(
                "Sample validation", CheckStatus.SKIP,
                "Could not read sample_list", group=group))
            return

        # Required columns
        required_cols = ['sample_name_unique', 'fastq_mate1_filename',
                         'fastq_mate2_filename', 'P5_index', 'P7_index']
        missing_cols = [c for c in required_cols if c not in sdf.columns]
        if missing_cols:
            report.add(CheckResult(
                "Required columns (sample_list)", CheckStatus.FAIL,
                f"Missing: {', '.join(missing_cols)}", group=group))
        else:
            report.add(CheckResult(
                "Required columns (sample_list)", CheckStatus.PASS, group=group))

        # Unique sample names
        if 'sample_name_unique' in sdf.columns:
            dups = sdf['sample_name_unique'].duplicated()
            if dups.any():
                dup_names = sdf.loc[dups, 'sample_name_unique'].tolist()[:5]
                report.add(CheckResult(
                    "Unique sample names", CheckStatus.FAIL,
                    f"{dups.sum()} duplicates: {', '.join(str(d) for d in dup_names)}",
                    group=group))
            else:
                report.add(CheckResult(
                    "Unique sample names", CheckStatus.PASS,
                    f"{len(sdf)} samples", group=group))

        # Duplicate barcodes within file pairs
        if 'barcode_name' in sdf.columns:
            file_pair_cols = [c for c in ['fastq_mate1_filename', 'fastq_mate2_filename']
                              if c in sdf.columns]
            if file_pair_cols:
                bc_col = 'barcode_name'
                grouped = sdf.groupby(file_pair_cols)[bc_col]
                dup_pairs = []
                for keys, grp in grouped:
                    if grp.duplicated().any():
                        dup_pairs.append(str(keys))
                if dup_pairs:
                    report.add(CheckResult(
                        "Duplicate barcodes", CheckStatus.FAIL,
                        f"Duplicate barcodes in file pairs: {', '.join(dup_pairs[:3])}",
                        group=group))
                else:
                    report.add(CheckResult(
                        "Duplicate barcodes", CheckStatus.PASS, group=group))

        # Index list checks
        # Format: columns 'type', 'id', 'sequence' with type values
        # like 'P5_index', 'P7_index', 'barcode'
        if index_list and Path(index_list).is_file():
            try:
                idf = pd.read_excel(index_list)
                idx_required = ['type', 'id', 'sequence']
                idx_missing = [c for c in idx_required if c not in idf.columns]
                if idx_missing:
                    report.add(CheckResult(
                        "Required columns (index_list)", CheckStatus.FAIL,
                        f"Missing: {', '.join(idx_missing)}", group=group))
                else:
                    report.add(CheckResult(
                        "Required columns (index_list)", CheckStatus.PASS,
                        group=group))

                # Check P5/P7/barcode IDs in sample_list match index_list
                if 'type' in idf.columns and 'id' in idf.columns:
                    index_ids = {}  # {type: set of ids}
                    for t in idf['type'].unique():
                        index_ids[t] = set(idf.loc[idf['type'] == t, 'id'].values)

                    unmatched = []
                    for col, idx_type in [('P5_index', 'P5_index'),
                                          ('P7_index', 'P7_index'),
                                          ('barcode', 'barcode')]:
                        if col in sdf.columns and idx_type in index_ids:
                            sample_vals = set(sdf[col].dropna().unique())
                            missing = sample_vals - index_ids[idx_type]
                            if missing:
                                unmatched.append(
                                    f"{col}: {', '.join(str(x) for x in sorted(missing)[:3])}")

                    if unmatched:
                        report.add(CheckResult(
                            "Index/barcode matching", CheckStatus.FAIL,
                            "IDs in sample_list not found in index_list:\n"
                            + "\n".join(f"  {u}" for u in unmatched),
                            group=group))
                    else:
                        report.add(CheckResult(
                            "Index/barcode matching", CheckStatus.PASS,
                            group=group))
            except Exception as e:
                report.add(CheckResult(
                    "Index list validation", CheckStatus.FAIL,
                    f"Cannot read: {e}", group=group))

        # Species match tRNA_database keys
        if 'species' in sdf.columns:
            tRNA_db = self.config.get('tRNA_database', {})
            if isinstance(tRNA_db, dict):
                sample_species = set(sdf['species'].dropna().unique())
                db_species = set(tRNA_db.keys())
                unmatched = sample_species - db_species
                if unmatched:
                    report.add(CheckResult(
                        "Species vs tRNA_database", CheckStatus.FAIL,
                        f"Species in sample_list with no database entry: "
                        f"{', '.join(sorted(unmatched))}",
                        group=group))
                else:
                    report.add(CheckResult(
                        "Species vs tRNA_database", CheckStatus.PASS,
                        group=group))

        # common_seqs single-species guard
        common_seqs = self.config.get('common_seqs')
        if common_seqs is not None and 'species' in sdf.columns:
            n_species = sdf['species'].nunique()
            if n_species > 1:
                report.add(CheckResult(
                    "common_seqs species guard", CheckStatus.WARN,
                    f"{n_species} species found but common_seqs is set — "
                    f"common_seqs is single-species only",
                    group=group))

        # abundance_control matches a sample_name
        if '7' in stages and self.config.get('run_abundance_analysis', False):
            ctrl = self.config.get('abundance_control')
            if ctrl and 'sample_name' in sdf.columns:
                if ctrl not in sdf['sample_name'].values:
                    report.add(CheckResult(
                        "abundance_control match", CheckStatus.FAIL,
                        f"'{ctrl}' not found in sample_name column",
                        group=group))
                else:
                    report.add(CheckResult(
                        "abundance_control match", CheckStatus.PASS,
                        group=group))

    def _preflight_tools(self, report, stages):
        """Check that required external tools are on PATH."""
        group = "Tools"

        # AdapterRemoval (stage 0a)
        if '0a' in stages:
            if shutil.which('AdapterRemoval'):
                report.add(CheckResult(
                    "AdapterRemoval", CheckStatus.PASS, group=group))
            else:
                report.add(CheckResult(
                    "AdapterRemoval", CheckStatus.FAIL,
                    "Not found on PATH. Install or load module.",
                    group=group))
        else:
            report.add(CheckResult(
                "AdapterRemoval", CheckStatus.SKIP,
                "stage 0a not selected", group=group))

        # swipe (stage 1)
        if '1' in stages:
            if shutil.which('swipe'):
                report.add(CheckResult("swipe", CheckStatus.PASS, group=group))
            else:
                report.add(CheckResult(
                    "swipe", CheckStatus.FAIL,
                    "Not found on PATH. Install or load module.",
                    group=group))
        else:
            report.add(CheckResult(
                "swipe", CheckStatus.SKIP,
                "stage 1 not selected", group=group))

    def log(self, message, level="INFO"):
        """Log message to file and console"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_msg = f"[{timestamp}] {level}: {message}"

        # Print to console
        print(log_msg)

        # Write to log file
        with open(self.log_file, 'a') as f:
            f.write(log_msg + '\n')

    def load_sample_info(self):
        """Load sample list and index information"""
        self.log("Loading sample information...")

        # Load sample list
        sample_list = self.config['sample_list']
        self.sample_df = pd.read_excel(sample_list)

        # Load index list
        index_list = self.config['index_list']
        index_df = pd.read_excel(index_list)

        # Add barcode and index sequences to sample_df
        self.sample_df = index_to_sample_df(self.sample_df, index_df)

        # Create input file DataFrame
        self.inp_file_df = self.sample_df[[
            'fastq_mate1_filename', 'fastq_mate2_filename',
            'P5_index', 'P7_index', 'P5_index_seq', 'P7_index_seq'
        ]].drop_duplicates().reset_index(drop=True)

        self.log(f"Loaded {len(self.sample_df)} samples")
        self.log(f"Input files: {len(self.inp_file_df)}")

    def _ensure_dir(self, *parts) -> Path:
        """Create and return a subdirectory under project_dir."""
        p = self.project_dir.joinpath(*parts)
        p.mkdir(parents=True, exist_ok=True)
        return p

    def setup_directories(self):
        """Setup directory structure"""
        self.log("Setting up directories...")

        self.dir_dict = {
            'NBdir': str(self.project_dir),
            'data_dir': 'data',
            'seq_dir': self.config.get('seq_dir', 'raw_fastq'),
            'AdapterRemoval_dir': 'AdapterRemoval',
            'BC_dir': 'BC_split',
            'UMI_dir': 'UMI_trimmed',
            'align_dir': 'SWalign',
            'stats_dir': 'stats_collection',
        }

        # Create data directory and preprocessing subdirectories
        (self.project_dir / 'data').mkdir(exist_ok=True)
        for dir_name in ['AdapterRemoval', 'BC_split', 'UMI_trimmed',
                         'SWalign', 'stats_collection']:
            (self.project_dir / 'data' / dir_name).mkdir(exist_ok=True)

        # Warn about legacy directories from previous versions
        for old in ['charge_analysis', 'fragment_analysis', 'modification_analysis',
                     'abundance_analysis', 'parquet_data']:
            if (self.project_dir / old).exists():
                self.log(f"  NOTE: Legacy directory '{old}/' found — new outputs go to results/", level="WARN")

    def _get_working_sample_df(self):
        """Return sample_df subset if --sample-index is set, else full."""
        if self.sample_index is not None:
            return self.sample_df.iloc[[self.sample_index]].copy()
        return self.sample_df

    def stage_0a_adapter_removal(self):
        """Stage 0a: Adapter Removal & Merging"""
        self.log("=" * 60)
        self.log("Stage 0a: Adapter Removal & Merging")
        self.log("=" * 60)

        MIN_READ_LEN = self.config.get('min_read_len', 39)

        # Run AdapterRemoval
        AR_obj = AR_merge(
            self.dir_dict, self.inp_file_df, MIN_READ_LEN,
            AR_threads=self.threads_per_job,
            overwrite_dir=self.config.get('overwrite', True)
        )
        self.inp_file_df = AR_obj.run_parallel(
            n_jobs=self.n_jobs,
            overwrite=self.config.get('overwrite', True)
        )

        # Log statistics
        for _, row in self.inp_file_df.iterrows():
            self.log(f"  {row['fastq_mate1_filename']}: "
                    f"{row['percent_successfully_merged']:.1f}% merged")

        # --- Merge QC probe: check UMI + barcode detectability ---
        try:
            barcodes = list(self.sample_df['barcode_seq'].unique())
            umi_mode = self.config.get('umi_trim_mode', 'anchored')
            probe_kwargs = {'umi_mode': umi_mode}

            if umi_mode == 'anchored':
                adapter_seqs_path = self.config.get('adapter_sequences', None)
                if adapter_seqs_path is None:
                    adapter_seqs_path = str(self.project_dir / 'utils' / 'adapter_sequences.yaml')
                with open(adapter_seqs_path) as f:
                    adapter_seqs = yaml.safe_load(f)
                anchor_name = self.config.get('umi_anchor')
                if anchor_name and anchor_name in adapter_seqs:
                    probe_kwargs['anchor_seq'] = adapter_seqs[anchor_name]
                    probe_kwargs['max_stagger'] = self.config.get('umi_max_stagger', 3)
                    probe_kwargs['anchor_max_dist'] = self.config.get('umi_anchor_max_dist', 1)

            probe_results = AR_obj.qc_probe(barcodes, **probe_kwargs)

            self.log(f"  Merge QC probe ({self.config.get('qc_probe_reads', 10000)} reads per file pair):")
            total_umi = total_bc = total_both = total_n = 0
            for basename, stats in probe_results.items():
                stagger_str = ""
                if 'stagger_counts' in stats:
                    counts = stats['stagger_counts']
                    stagger_str = f"  stagger=[{', '.join(str(counts.get(i, 0)) for i in range(max(counts.keys()) + 1))}]"
                self.log(f"    {basename}: UMI={stats['pct_umi']}% "
                         f"BC={stats['pct_bc']}% both={stats['pct_both']}%{stagger_str}")
                total_umi += stats['pct_umi'] * stats['n_probed']
                total_bc += stats['pct_bc'] * stats['n_probed']
                total_both += stats['pct_both'] * stats['n_probed']
                total_n += stats['n_probed']

            if total_n > 0:
                avg_umi = total_umi / total_n
                avg_bc = total_bc / total_n
                avg_both = total_both / total_n
                self.log(f"  Average: UMI={avg_umi:.1f}% BC={avg_bc:.1f}% both={avg_both:.1f}%")
                if avg_umi < 50:
                    self.log("  WARNING: Low UMI detection rate (<50%). "
                             "Check umi_trim_mode / umi_anchor config.", level="WARN")
                if avg_bc < 50:
                    self.log("  WARNING: Low barcode detection rate (<50%). "
                             "Check barcode sequences in sample list.", level="WARN")
        except Exception as e:
            self.log(f"  Merge QC probe skipped: {e}", level="WARN")

        self.status['stages_completed'].append('0a')

    def stage_0b_barcode_split(self):
        """Stage 0b: Barcode Splitting"""
        self.log("=" * 60)
        self.log("Stage 0b: Barcode Splitting")
        self.log("=" * 60)

        # Run barcode splitting
        BCsplit_obj = BC_split(
            self.dir_dict, self.sample_df, self.inp_file_df,
            overwrite_dir=self.config.get('overwrite', True)
        )
        self.sample_df, self.inp_file_df = BCsplit_obj.run_parallel(
            n_jobs=self.n_jobs
        )

        # Log statistics
        for _, row in self.inp_file_df.iterrows():
            if row['percent_BC-mapped'] < 80:
                self.log(f"  WARNING: Low barcode mapping: "
                        f"{row['percent_BC-mapped']:.1f}%", level="WARN")
            else:
                self.log(f"  Barcode mapping: {row['percent_BC-mapped']:.1f}%")

        self.status['stages_completed'].append('0b')

    def stage_0c_umi_trim(self):
        """Stage 0c: UMI Trimming"""
        self.log("=" * 60)
        self.log("Stage 0c: UMI Trimming")
        self.log("=" * 60)

        # Run UMI trimming
        downsample_percentile = self.config.get('downsample_percentile', None)
        downsample_absolute = self.config.get('downsample_absolute', None)

        # UMI trim mode: 'anchored' (default) or 'pyrimidine' (legacy)
        umi_mode = self.config.get('umi_trim_mode', 'anchored')
        anchor_seq = None
        if umi_mode == 'anchored':
            adapter_seqs_path = self.config.get('adapter_sequences', None)
            if adapter_seqs_path is None:
                adapter_seqs_path = str(self.project_dir / 'utils' / 'adapter_sequences.yaml')
            with open(adapter_seqs_path) as f:
                adapter_seqs = yaml.safe_load(f)
            anchor_name = self.config.get('umi_anchor')
            if anchor_name is None:
                raise ValueError("Config field 'umi_anchor' is required when umi_trim_mode='anchored'")
            if anchor_name not in adapter_seqs:
                raise ValueError(f"Anchor '{anchor_name}' not found in {adapter_seqs_path}. "
                                 f"Available: {list(adapter_seqs.keys())}")
            anchor_seq = adapter_seqs[anchor_name]
            self.log(f"  UMI mode: anchored (anchor={anchor_name}: {anchor_seq})")
        else:
            self.log(f"  UMI mode: pyrimidine (legacy)")

        working_df = self._get_working_sample_df()
        UMItrim_obj = UMI_trim(
            self.dir_dict, working_df,
            mode=umi_mode,
            anchor_seq=anchor_seq,
            max_stagger=self.config.get('umi_max_stagger', 3),
            anchor_max_dist=self.config.get('umi_anchor_max_dist', 1),
            overwrite_dir=self.config.get('overwrite', True) if self.sample_index is None else False,
            downsample_percentile=downsample_percentile,
            downsample_absolute=downsample_absolute
        )
        result_df = UMItrim_obj.run_parallel(n_jobs=self.n_jobs)
        if self.sample_index is not None:
            # Merge single-sample results back into full sample_df
            merge_cols = [c for c in result_df.columns if c not in self.sample_df.columns or c == 'sample_name_unique']
            extra_cols = [c for c in result_df.columns if c not in self.sample_df.columns]
            if extra_cols:
                self.sample_df = self.sample_df.merge(
                    result_df[['sample_name_unique'] + extra_cols],
                    on='sample_name_unique', how='left')
        else:
            self.sample_df = result_df

        # Log statistics
        avg_umi_pct = self.sample_df['percent_UMI_obs-vs-exp'].mean()
        self.log(f"  Average UMI obs/exp: {avg_umi_pct:.1f}%")

        if avg_umi_pct < 80:
            self.log(f"  WARNING: Low UMI diversity suggests library prep bottleneck",
                    level="WARN")

        self.status['stages_completed'].append('0c')

    def stage_1_alignment(self):
        """Stage 1: SWIPE Alignment"""
        self.log("=" * 60)
        self.log("Stage 1: SWIPE Alignment")
        self.log("=" * 60)

        # Load tRNA database info
        tRNA_database = self.config['tRNA_database']
        self.tRNA_data = read_tRNAdb_info(tRNA_database)

        # Setup alignment parameters
        SWIPE_score_mat = self.config['SWIPE_score_mat']
        common_seqs = self.config.get('common_seqs', None)
        MIN_SCORE_ALIGN = self.config.get('min_score_align', 15)

        # Run alignment
        working_df = self._get_working_sample_df()
        align_obj = SWIPE_align(
            self.dir_dict, tRNA_database, working_df,
            SWIPE_score_mat,
            gap_penalty=self.config.get('gap_penalty', 6),
            extension_penalty=self.config.get('extension_penalty', 3),
            min_score_align=MIN_SCORE_ALIGN,
            common_seqs=common_seqs,
            overwrite_dir=self.config.get('overwrite', True) if self.sample_index is None else False,
            SWIPE_threads=self.threads_per_job,
        )
        result_df = align_obj.run_parallel(n_jobs=self.n_jobs)
        if self.sample_index is not None:
            merge_cols = [c for c in result_df.columns if c not in self.sample_df.columns]
            if merge_cols:
                self.sample_df = self.sample_df.merge(
                    result_df[['sample_name_unique'] + merge_cols],
                    on='sample_name_unique', how='left')
        else:
            self.sample_df = result_df

        # Log statistics
        avg_mapping = self.sample_df['Mapping_percent'].mean()
        self.log(f"  Average mapping rate: {avg_mapping:.1f}%")

        if avg_mapping < 80:
            self.log(f"  WARNING: Low mapping rate", level="WARN")

        self.status['stages_completed'].append('1')

    def stage_2_stats_collection(self):
        """Stage 2: Stats Collection"""
        self.log("=" * 60)
        self.log("Stage 2: Stats Collection")
        self.log("=" * 60)

        common_seqs = self.config.get('common_seqs', None)

        align_dir = self.project_dir / 'data' / self.dir_dict['align_dir']
        common_obs_files = list(align_dir.glob('*_common-seq-obs.json*'))

        if common_seqs is not None and not common_obs_files:
            # common_seqs configured but no result files from alignment
            self.log("  WARNING: common_seqs configured but no common-seq-obs files "
                    "found in SWalign/. Setting common_seqs=None.", level="WARN")
            common_seqs = None
        elif common_seqs is None and common_obs_files:
            # Stale common-seq-obs files from a previous run with common_seqs
            self.log(f"  WARNING: Removing {len(common_obs_files)} stale common-seq-obs "
                     "files from SWalign/ (common_seqs is now null).", level="WARN")
            for f in common_obs_files:
                f.unlink()
            common_obs_files = []

        # Run stats collection
        stats_obj = STATS_collection(
            self.dir_dict, self.tRNA_data, self.sample_df,
            common_seqs=common_seqs,
            overwrite_dir=self.config.get('overwrite', True)
        )
        self.stats_df = stats_obj.run_parallel(n_jobs=self.n_jobs)

        # Copy fragment counts to sample_df
        self.fragment_counts = getattr(stats_obj, 'fragment_counts', {})
        if self.fragment_counts:
            for col in ['N_full_length', 'N_rt_dropoff', 'N_5p_fragment', 'N_degraded', 'N_total_aligned']:
                self.sample_df[col] = self.sample_df['sample_name_unique'].map(
                    lambda sn, c=col: self.fragment_counts.get(sn, {}).get(c, 0)
                )

            # Recalculate Mapping_percent from fragment counts:
            # N_total_aligned / N_after_trim * 100
            if 'N_after_trim' in self.sample_df.columns:
                self.sample_df['Mapping_percent'] = np.where(
                    self.sample_df['N_after_trim'] > 0,
                    self.sample_df['N_total_aligned'] / self.sample_df['N_after_trim'] * 100,
                    0.0
                )

        self.log(f"  Collected stats for {len(self.stats_df)} entries")

        self.status['stages_completed'].append('2')

    def stage_3_charge_quantification(self):
        """Stage 3: Charge Quantification"""
        self.log("=" * 60)
        self.log("Stage 3: Charge Quantification")
        self.log("=" * 60)

        if not CHARGE_AVAILABLE:
            self.log("WARNING: ChargeQuantifier not available. Skipping stage 3.", level="WARN")
            self.log("  Install with: pip install -e . (from repo root)", level="WARN")
            return

        # Get stats file
        stats_file = self.project_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
        if not stats_file.exists():
            self.log(f"ERROR: Stats file not found: {stats_file}", level="ERROR")
            return

        # Get charge quantification settings
        charge_count = self.config.get('charge_count', 'count')
        charge_levels = self.config.get('charge_levels', ['transcript', 'codon', 'aa'])
        include_synthetic = self.config.get('include_synthetic', False)
        include_mito = self.config.get('include_mito', True)

        self.log(f"  Input: {stats_file.name}")
        self.log(f"  Charge count column: {charge_count}")
        self.log(f"  Levels: {', '.join(charge_levels)}")

        try:
            # tRNA source classification prefixes
            source_prefixes = self.config.get(
                'tRNA_source_prefixes', {'Synthetic_': 'synthetic'})

            # Initialize quantifier
            quantifier = ChargeQuantifier(
                stats_csv=str(stats_file),
                charge_count=charge_count,
                RPM_count=charge_count,
                source_prefixes=source_prefixes,
            )

            self.log(f"  Loaded {len(quantifier.stats_df)} alignment records")

            # Create output directory
            charge_dir = self._ensure_dir('results', 'charge')

            # Quantify and export for each level
            self.charge_results = {}
            for level in charge_levels:
                self.log(f"  Quantifying at {level} level...")

                df = quantifier.quantify_all(
                    level=level,
                    include_synthetic=include_synthetic,
                    include_mito=include_mito
                )

                # Save to CSV
                output_file = charge_dir / f'charge_df_{level}.csv'
                df.to_csv(output_file, index=False)
                self.charge_results[level] = df

                self.log(f"    Saved: {output_file.name} ({len(df)} entries)")

            # Generate summary statistics
            self.log("  Generating summary statistics...")
            summary_list = []
            for level in charge_levels:
                summary = quantifier.get_summary_statistics(level=level)
                summary['level'] = level
                summary_list.append(summary)

            if summary_list:
                summary_df = pd.concat(summary_list, ignore_index=True)
                summary_file = charge_dir / 'charge_summary.csv'
                summary_df.to_csv(summary_file, index=False)
                self.charge_summary_df = summary_df
                self.log(f"    Saved: {summary_file.name}")

            self.status['stages_completed'].append('3')
            self.log("  Charge quantification complete!")

            # Generate charge report
            if CHARGE_REPORT_AVAILABLE:
                try:
                    sample_df = getattr(self, 'sample_df', None)
                    charge_tr = self.charge_results.get('transcript')
                    charge_aa = self.charge_results.get('aa')
                    if charge_tr is not None and charge_aa is not None and sample_df is not None:
                        source_prefixes = self.config.get('tRNA_source_prefixes',
                                                          {'Synthetic_': 'synthetic'})
                        report_gen = ChargeReportGenerator(
                            charge_df_transcript=charge_tr,
                            charge_df_aa=charge_aa,
                            charge_summary=self.charge_summary_df,
                            sample_df=sample_df,
                            source_prefixes=source_prefixes,
                        )
                        qc_dir = self.project_dir / 'qc_reports'
                        qc_dir.mkdir(parents=True, exist_ok=True)
                        report_path = qc_dir / 'charge_report.html'
                        report_gen.generate_html_report(report_path)
                        self.log(f"  Saved: {report_path.name}")
                except Exception as report_err:
                    self.log(f"  WARNING: Could not generate charge report: "
                             f"{report_err}", level="WARN")

        except Exception as e:
            self.log(f"ERROR in charge quantification: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def _run_fragment_analysis(self):
        """Run fragment analysis as part of stage 3.

        Reads per-sample stats CSVs and classifies reads into fragment
        types, profiles RT drop-off positions, and computes length
        distributions.
        """
        if not FRAGMENTS_AVAILABLE:
            self.log("  Fragment analysis skipped (module not available)", level="WARN")
            return

        stats_dir = self.project_dir / 'data' / 'stats_collection'
        if not stats_dir.exists():
            self.log("  Fragment analysis skipped (no stats_collection dir)", level="WARN")
            return

        sample_names = self.sample_df['sample_name_unique'].tolist()
        self.log(f"  Running fragment analysis ({len(sample_names)} samples)...")

        try:
            analyser = FragmentAnalyser(
                stats_dir=stats_dir,
                sample_names=sample_names,
                min_reads=self.config.get('fragment_min_reads', 10),
            )
            analyser.run()

            frag_dir = self._ensure_dir('results', 'fragments')
            write_csv = self.config.get('fragment_write_csv', False)
            analyser.export(frag_dir, write_csv=write_csv)

            summary = analyser._summary
            if summary is not None and not summary.empty:
                n_samples = len(summary)
                fl = summary['pct_full_length']
                rt = summary['pct_rt_dropoff']
                self.log(f"  Fragment analysis complete ({n_samples} samples):")
                self.log(f"    Full-length:  {fl.min():.1f} – {fl.max():.1f}%  (mean {fl.mean():.1f}%)")
                self.log(f"    RT drop-off:  {rt.min():.1f} – {rt.max():.1f}%  (mean {rt.mean():.1f}%)")
                self.log(f"    See results/fragments/ and qc_reports/fragment_report.html for details")
            else:
                self.log("  Fragment analysis complete!")

            self._fragment_analysis_ran = True

            # Generate fragment report
            if FRAGMENT_REPORT_AVAILABLE:
                try:
                    sample_df = getattr(self, 'sample_df', None)
                    if sample_df is not None:
                        source_prefixes = self.config.get('tRNA_source_prefixes',
                                                          {'Synthetic_': 'synthetic'})
                        coverage_df = analyser._coverage if analyser._coverage is not None else None
                        report_gen = FragmentReportGenerator(
                            fragment_counts_df=analyser._fragment_counts if analyser._fragment_counts is not None else pd.DataFrame(),
                            rt_dropoff_df=analyser._rt_dropoff if analyser._rt_dropoff is not None else pd.DataFrame(),
                            fragment_lengths_df=analyser._fragment_lengths if analyser._fragment_lengths is not None else pd.DataFrame(),
                            fragment_summary_df=summary if summary is not None else pd.DataFrame(),
                            sample_df=sample_df,
                            coverage_df=coverage_df,
                            source_prefixes=source_prefixes,
                        )
                        qc_dir = self.project_dir / 'qc_reports'
                        qc_dir.mkdir(parents=True, exist_ok=True)
                        report_path = qc_dir / 'fragment_report.html'
                        report_gen.generate_html_report(report_path)
                        self.log(f"  Saved: {report_path.name}")
                except Exception as report_err:
                    self.log(f"  WARNING: Could not generate fragment report: "
                             f"{report_err}", level="WARN")

        except Exception as e:
            self.log(f"  ERROR in fragment analysis: {e}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def stage_4_parquet_storage(self):
        """Stage 4: Parquet Storage (optional)"""
        self.log("=" * 60)
        self.log("Stage 4: Parquet Storage")
        self.log("=" * 60)

        if not STORAGE_AVAILABLE:
            self.log("WARNING: tRNAseqDataStore not available. Skipping stage 4.", level="WARN")
            self.log("  Install with: pip install pyarrow", level="WARN")
            return

        try:
            # Initialize data store
            parquet_dir = self._ensure_dir('results', 'parquet')
            store = tRNAseqDataStore(str(parquet_dir))

            compression = self.config.get('parquet_compression', 'snappy')
            self.log(f"  Compression: {compression}")

            # Convert stats CSV to Parquet
            stats_file = self.project_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
            if stats_file.exists():
                self.log(f"  Converting stats to Parquet...")
                stats_df = pd.read_csv(stats_file)
                parquet_file = store.save_stats(
                    stats_df,
                    filename='ALL_stats_aggregate',
                    format='parquet',
                    compression=compression
                )

                # Calculate compression ratio
                csv_size = stats_file.stat().st_size / (1024 * 1024)
                parquet_size = parquet_file.stat().st_size / (1024 * 1024)
                ratio = csv_size / parquet_size if parquet_size > 0 else 0

                self.log(f"    CSV: {csv_size:.1f} MB")
                self.log(f"    Parquet: {parquet_size:.1f} MB")
                self.log(f"    Compression ratio: {ratio:.1f}x")

            # Convert charge CSVs to Parquet if they exist
            charge_dir = self.project_dir / 'results' / 'charge'
            if charge_dir.exists():
                charge_parquet_dir = parquet_dir / 'charge'
                charge_parquet_dir.mkdir(parents=True, exist_ok=True)
                for csv_file in charge_dir.glob('charge_df_*.csv'):
                    self.log(f"  Converting {csv_file.name} to Parquet...")
                    df = pd.read_csv(csv_file)
                    out_path = charge_parquet_dir / csv_file.with_suffix('.parquet').name
                    df.to_parquet(out_path, compression=compression)
                    self.log(f"    Saved: {out_path.name}")

            # Copy fragment parquet files into results/parquet/fragments/
            frag_dir = self.project_dir / 'results' / 'fragments'
            if frag_dir.exists():
                frag_parquet_dir = parquet_dir / 'fragments'
                frag_parquet_dir.mkdir(parents=True, exist_ok=True)
                for pq_file in frag_dir.glob('*.parquet'):
                    import shutil
                    shutil.copy(pq_file, frag_parquet_dir / pq_file.name)
                    self.log(f"  Copied fragment: {pq_file.name}")

            self.status['stages_completed'].append('4')
            self.log("  Parquet conversion complete!")

        except Exception as e:
            self.log(f"ERROR in Parquet storage: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def stage_5_qc_report(self):
        """Stage 5: QC Summary Report

        Generates QC_summary.csv and QC_report.html from data already
        computed by earlier stages. No new computation — just formatting.
        """
        self.log("=" * 60)
        self.log("Stage 5: QC Summary Report")
        self.log("=" * 60)

        if not QC_AVAILABLE:
            self.log("WARNING: QCReportGenerator not available. Skipping stage 5.", level="WARN")
            return

        try:
            # Load stats_df if not already in memory
            stats_df = getattr(self, 'stats_df', None)
            if stats_df is None:
                stats_file = (self.project_dir / 'data' / 'stats_collection'
                              / 'ALL_stats_aggregate.csv')
                if stats_file.exists():
                    stats_df = pd.read_csv(stats_file)

            # Load charge summary if not already in memory
            charge_summary = getattr(self, 'charge_summary_df', None)
            if charge_summary is None:
                cs_file = self.project_dir / 'results' / 'charge' / 'charge_summary.csv'
                if cs_file.exists():
                    charge_summary = pd.read_csv(cs_file)

            # Ensure we have sample_df and inp_file_df
            sample_df = getattr(self, 'sample_df', None)
            inp_file_df = getattr(self, 'inp_file_df', pd.DataFrame())

            if sample_df is None:
                self.log("ERROR: sample_df not available for QC report", level="ERROR")
                return

            qc = QCReportGenerator(
                project_dir=self.project_dir,
                sample_df=sample_df,
                inp_file_df=inp_file_df,
                charge_summary_df=charge_summary,
                stats_df=stats_df,
                bc_dir=self.project_dir / 'data' / 'BC_split',
            )

            # QC outputs go into qc_reports/
            qc_dir = self.project_dir / 'qc_reports'
            qc_dir.mkdir(exist_ok=True)

            summary_path = qc_dir / 'QC_summary.csv'
            qc_summary = qc.save_summary_csv(summary_path)
            self.log(f"  Saved: qc_reports/{summary_path.name} ({len(qc_summary)} samples)")

            report_path = qc_dir / 'QC_report.html'
            qc.generate_html_report(report_path)
            self.log(f"  Saved: qc_reports/{report_path.name}")

            self.status['stages_completed'].append('5')
            self.log("  QC report generation complete!")

        except Exception as e:
            self.log(f"ERROR in QC report: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def stage_6_modification_analysis(self):
        """Stage 6: Modification Analysis (optional)

        Six phases:
        1. PSCM extraction (per-sample positional count matrices)
        2. Background error-rate estimation (synthetic spike-ins or empirical)
        3. Per-sample modification calling (RT signatures + binomial test)
        4. Replicate aggregation (Fisher combined p-value + consensus sieve)
        5. Summary CSV (one row per sample)
        6. QC report (interactive HTML dashboard)
        """
        self.log("=" * 60)
        self.log("Stage 6: Modification Analysis")
        self.log("=" * 60)

        if not MODIFICATIONS_AVAILABLE:
            self.log("WARNING: Modification analysis modules not available. "
                     "Skipping stage 6.", level="WARN")
            return

        try:
            # ---- Config ----
            tRNA_database = self.config['tRNA_database']
            ref_fasta = list(tRNA_database.values())[0]

            organism = self.config.get('organism', 'Escherichia coli')
            min_coverage = self.config.get('modification_min_coverage', 50)
            discover_novel = self.config.get('discover_novel_modifications', False)
            use_api = not self.config.get('no_modomics', True)
            min_replicates = self.config.get('modification_min_replicates', 3)
            mod_alpha = self.config.get('modification_alpha', 0.01)
            synthetic_prefixes = tuple(
                self.config.get('synthetic_tRNA_prefixes', ['Synthetic_'])
            )

            json_dir = self.project_dir / 'data' / self.dir_dict['align_dir']
            output_dir = self._ensure_dir('results', 'modifications')

            sample_names = self.sample_df['sample_name_unique'].tolist()

            self.log(f"  Reference: {ref_fasta}")
            self.log(f"  Organism: {organism}")
            self.log(f"  Samples: {len(sample_names)}")
            self.log(f"  Novel discovery: {discover_novel}")

            # ==== Phase 1: PSCM extraction ====
            self.log("  Phase 1/7: Extracting PSCMs...")
            extractor = PositionalExtractor(ref_fasta)
            all_pscm = extractor.run_parallel(
                json_dir, sample_names, n_jobs=self.n_jobs
            )

            # ==== Phase 2: Background estimation ====
            self.log("  Phase 2/7: Estimating background error rate...")
            pooled = _pool_pscm_dicts(all_pscm)
            bg_rate, bg_source = estimate_background_error_rate(
                pooled, extractor.ref_dict,
                synthetic_prefixes=synthetic_prefixes,
                min_coverage=min_coverage,
            )
            self.log(f"  Background error rate: {bg_rate:.5f} (source: {bg_source})")

            # MODOMICS
            annotator = MODOMICSAnnotator(organism)
            mods_df = annotator.get_modifications(use_api=use_api)
            self.log(f"  Loaded {len(mods_df)} known modification entries")

            # ==== Phase 3: Per-sample calling ====
            self.log("  Phase 3/7: Calling modifications per sample...")
            analyzer = RTSignatureAnalyzer(
                min_coverage=min_coverage, verbose=False
            )
            analyzer.load_reference(ref_fasta)
            caller = ModificationCaller(
                organism=organism,
                background_error_rate=bg_rate,
                alpha=mod_alpha,
            )

            per_sample_calls = {}
            sample_call_counts = []

            for sample_name, pscm_dict in all_pscm.items():
                pscm_dfs = analyzer.load_pscm_from_positional(pscm_dict)

                rt_profile = extractor.compute_rt_profile(pscm_dict)
                mm_profile = extractor.compute_mismatch_profile(pscm_dict)

                sample_dir = output_dir / sample_name
                sample_dir.mkdir(exist_ok=True)

                _save_df(rt_profile, sample_dir / 'rt_profile')
                _save_df(mm_profile, sample_dir / 'mismatch_profile')

                # Detect anticodon positions (fallback for CSV-only data)
                ac_positions = extractor._autodetect_anticodon_positions(
                    list(pscm_dfs.keys())
                )

                # Call modifications for each tRNA in this sample
                sample_calls = []
                for trna_name, pscm_df in pscm_dfs.items():
                    rt_counts = pscm_dict[trna_name][:, 7]
                    analysis = analyzer.analyze_trna_with_actual_stops(
                        trna_name, pscm_df, rt_stop_counts=rt_counts,
                    )
                    ref_info = analyzer.reference_sequences.get(trna_name, {})
                    ref_seq = ref_info.get('seq')

                    # Get known mods: alignment-based (preferred) or heuristic fallback
                    known_mods_for_trna = None
                    if ref_seq:
                        ac_pos = ac_positions.get(trna_name)
                        ac_start = ac_pos[0] if ac_pos is not None else None
                        known_mods_for_trna = annotator.get_known_mods_linear(
                            trna_name, ref_seq,
                            anticodon_linear_start=ac_start,
                        )

                    calls_df = caller.call_all(
                        trna_name,
                        analysis['signatures'],
                        pscm_df,
                        ref_seq,
                        discover_novel=discover_novel,
                        min_coverage=min_coverage,
                        known_mods_df=known_mods_for_trna,
                    )
                    if not calls_df.empty:
                        sample_calls.append(calls_df)

                if sample_calls:
                    combined = pd.concat(sample_calls, ignore_index=True)
                    per_sample_calls[sample_name] = combined
                    _save_df(combined, sample_dir / 'modification_calls')
                    sample_call_counts.append(len(combined))
                else:
                    per_sample_calls[sample_name] = pd.DataFrame()
                    sample_call_counts.append(0)

            # Log modification summary
            n_mod_samples = len(sample_call_counts)
            total_calls = sum(sample_call_counts)
            unique_trnas = len(set().union(*(
                set(df['trna_name']) for df in per_sample_calls.values()
                if not df.empty and 'trna_name' in df.columns
            ))) if per_sample_calls else 0
            if n_mod_samples > 0 and total_calls > 0:
                counts = np.array(sample_call_counts)
                self.log(f"  Modification calling complete ({n_mod_samples} samples):")
                self.log(f"    Total: {total_calls:,} calls across {unique_trnas} unique tRNAs")
                self.log(f"    Per-sample: {counts.min()} – {counts.max()} calls  (mean {counts.mean():.1f})")
                self.log(f"    See qc_reports/modification_report.html for per-sample details")
            else:
                self.log(f"  Modification calling complete ({n_mod_samples} samples): 0 calls")

            # ==== Phase 4: Replicate aggregation ====
            self.log("  Phase 4/7: Aggregating across replicates...")
            replicate_groups = {}
            if 'sample_name' in self.sample_df.columns:
                for _, row in self.sample_df.iterrows():
                    snu = row['sample_name_unique']
                    sn = str(row['sample_name'])
                    replicate_groups.setdefault(sn, []).append(snu)

            aggregated_calls = pd.DataFrame()
            consensus_calls = pd.DataFrame()

            if replicate_groups:
                aggregator = ReplicateAggregator(
                    min_replicates=min_replicates, alpha=mod_alpha
                )
                aggregated_calls = aggregator.aggregate(
                    per_sample_calls, replicate_groups
                )
                if not aggregated_calls.empty:
                    _save_df(aggregated_calls,
                             output_dir / 'aggregated_modifications')
                    consensus_calls = aggregated_calls[
                        aggregated_calls['consensus_call']
                    ].copy()
                    _save_df(consensus_calls,
                             output_dir / 'consensus_modifications')
                    self.log(f"  Aggregated: {len(aggregated_calls)} sites, "
                             f"{len(consensus_calls)} consensus")
                else:
                    self.log("  No aggregated modification calls.")
            else:
                self.log("  No replicate groups found — skipping aggregation.")

            # ==== Phase 5: Summary CSV ====
            self.log("  Phase 5/7: Generating modification summary...")
            summary_rows = []
            for snu in sample_names:
                sc = per_sample_calls.get(snu, pd.DataFrame())
                n_total = len(sc)
                n_consensus = 0
                if not consensus_calls.empty and not sc.empty:
                    consensus_sites = set(
                        zip(consensus_calls.get('trna_name', []),
                            consensus_calls.get('position', []))
                    )
                    sample_sites = set(zip(sc['trna_name'], sc['position']))
                    n_consensus = len(sample_sites & consensus_sites)

                mean_fc = float(sc['fold_change'].mean()) if (
                    not sc.empty and 'fold_change' in sc.columns
                ) else np.nan

                summary_rows.append({
                    'sample_name_unique': snu,
                    'total_calls': n_total,
                    'consensus_calls': n_consensus,
                    'mean_fold_change': round(mean_fc, 3) if not np.isnan(mean_fc) else np.nan,
                    'background_error_rate': bg_rate,
                    'bg_source': bg_source,
                })

            summary_df = pd.DataFrame(summary_rows)
            summary_df.to_csv(output_dir / 'modification_summary.csv', index=False)
            self.log(f"  Saved modification_summary.csv ({len(summary_df)} rows)")

            # ==== Phase 6: Crosstalk analysis (SLAC) ====
            self.log("  Phase 6/7: Analyzing modification crosstalks (SLAC)...")
            crosstalk_df = pd.DataFrame()
            try:
                from trnaseq.modifications.crosstalk import CrosstalkAnalyzer

                # Build mod_positions dict: {trna_name: [linear_pos, ...]}
                mod_pos_dict = {}
                for trna_name, ref_info in analyzer.reference_sequences.items():
                    ref_seq = ref_info.get('seq')
                    if ref_seq:
                        known = annotator.get_known_mods_linear(trna_name, ref_seq)
                        if not known.empty and 'linear_position' in known.columns:
                            positions = sorted(known['linear_position'].astype(int).tolist())
                            if len(positions) >= 2:
                                mod_pos_dict[trna_name] = positions

                if mod_pos_dict:
                    ct_analyzer = CrosstalkAnalyzer(
                        min_coverage=min_coverage,
                        mismatch_threshold=0.05,
                    )
                    # Collect SWalign JSON paths per sample
                    swalign_dir = self.project_dir / 'data' / 'SWalign'
                    json_paths = {}
                    for snu in sample_names:
                        jp = swalign_dir / f'{snu}_SWalign.json.bz2'
                        if jp.exists():
                            json_paths[snu] = jp

                    if json_paths:
                        crosstalk_df = ct_analyzer.analyze_multiple_samples(
                            json_paths, mod_pos_dict,
                            n_jobs=self.n_jobs,
                        )
                        if not crosstalk_df.empty:
                            _save_df(crosstalk_df, output_dir / 'crosstalk_analysis')
                            n_sig = crosstalk_df['fdr_significant'].sum() if 'fdr_significant' in crosstalk_df.columns else 0
                            n_pairs = crosstalk_df[['trna_name', 'pos_a', 'pos_b']].drop_duplicates().shape[0]
                            self.log(f"  Crosstalk: {n_pairs} position pairs tested, "
                                     f"{n_sig} significant (FDR < 0.05)")
                        else:
                            self.log("  No crosstalks detected (insufficient coverage or signal).")
                    else:
                        self.log("  No SWalign JSON files found — skipping crosstalk.")
                else:
                    self.log("  No tRNAs with ≥2 known modification positions — skipping crosstalk.")
            except ImportError:
                self.log("  Crosstalk module not available — skipping.")
            except Exception as ct_err:
                self.log(f"  WARNING: Crosstalk analysis failed: {ct_err}", level="WARN")

            # ==== Phase 7: QC report ====
            self.log("  Phase 7/7: Generating modification QC report...")
            if MOD_REPORT_AVAILABLE:
                try:
                    source_prefixes = self.config.get('tRNA_source_prefixes',
                                                      {'Synthetic_': 'synthetic'})
                    report_gen = ModificationReportGenerator(
                        per_sample_calls=per_sample_calls,
                        aggregated_calls=aggregated_calls if not aggregated_calls.empty else None,
                        consensus_calls=consensus_calls if not consensus_calls.empty else None,
                        replicate_groups=replicate_groups if replicate_groups else None,
                        ref_dict=extractor.ref_dict,
                        summary_df=summary_df,
                        source_prefixes=source_prefixes,
                    )
                    qc_dir = self.project_dir / 'qc_reports'
                    qc_dir.mkdir(parents=True, exist_ok=True)
                    report_path = report_gen.generate_html_report(
                        qc_dir / 'modification_report.html'
                    )
                    self.log(f"  Saved modification report: {report_path}")
                except Exception as report_err:
                    self.log(f"  WARNING: Could not generate modification report: "
                             f"{report_err}", level="WARN")
            else:
                self.log("  Modification report generator not available — skipping.")

            self.status['stages_completed'].append('6')
            self.log("  Modification analysis complete!")

        except Exception as e:
            self.log(f"ERROR in modification analysis: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def stage_7_abundance_analysis(self):
        """Stage 7: Differential Abundance Analysis (optional)

        Runs DESeq2 via pyDESeq2 to identify differentially abundant tRNAs
        between conditions. Generates results CSVs and an interactive HTML
        dashboard.
        """
        self.log("=" * 60)
        self.log("Stage 7: Differential Abundance Analysis")
        self.log("=" * 60)

        if not ABUNDANCE_AVAILABLE:
            self.log("WARNING: DifferentialAbundance not available. "
                     "Skipping stage 7.", level="WARN")
            self.log("  Install with: pip install pydeseq2", level="WARN")
            return

        try:
            stats_file = (self.project_dir / 'data' / 'stats_collection'
                          / 'ALL_stats_aggregate.csv')
            if not stats_file.exists():
                self.log(f"ERROR: Stats file not found: {stats_file}",
                         level="ERROR")
                return

            sample_df = getattr(self, 'sample_df', None)
            if sample_df is None:
                self.log("ERROR: sample_df not available", level="ERROR")
                return

            level = self.config.get('abundance_level', 'aa')
            control = self.config.get('abundance_control', None)

            self.log(f"  Level: {level}")
            self.log(f"  Control group: {control or '(auto-detect)'}")

            da = DifferentialAbundance(
                stats_csv=str(stats_file),
                sample_df=sample_df,
                level=level,
                control_group=control,
            )
            self.log(f"  Count matrix: {da.count_matrix.shape[0]} samples x "
                     f"{da.count_matrix.shape[1]} features")
            self.log(f"  Control: {da.control_group}")

            results = da.run_deseq2()
            self.log(f"  DESeq2 results: {len(results)} comparisons")

            # Export results
            abundance_dir = self._ensure_dir('results', 'abundance')
            da.export_results(abundance_dir)
            self.log(f"  Saved results to results/abundance/")

            # Generate report
            if ABUNDANCE_REPORT_AVAILABLE and not results.empty:
                try:
                    report_gen = AbundanceReportGenerator(
                        results_df=results,
                        count_matrix=da.count_matrix,
                        control_group=da.control_group,
                        level=level,
                        condition_map=da.condition_map,
                    )
                    qc_dir = self.project_dir / 'qc_reports'
                    qc_dir.mkdir(parents=True, exist_ok=True)
                    report_path = qc_dir / 'abundance_report.html'
                    report_gen.generate_html_report(report_path)
                    self.log(f"  Saved: abundance_report.html")
                except Exception as report_err:
                    self.log(f"  WARNING: Could not generate abundance report: "
                             f"{report_err}", level="WARN")

            self.status['stages_completed'].append('7')
            self.log("  Abundance analysis complete!")

        except Exception as e:
            self.log(f"ERROR in abundance analysis: {str(e)}", level="ERROR")
            import traceback
            self.log(traceback.format_exc(), level="ERROR")

    def _file_size_mb(self, path):
        """Return file size in MB, or None if file doesn't exist."""
        p = Path(path)
        if p.exists():
            return round(p.stat().st_size / (1024 * 1024), 3)
        return None

    def _collect_file_metrics(self):
        """Collect per-sample file sizes from pipeline output directories."""
        data_dir = self.project_dir / 'data'
        samples = self.sample_df['sample_name_unique'].tolist()

        for sample in samples:
            self.file_metrics[sample] = {
                'bc_split_fastq_MB': self._file_size_mb(
                    data_dir / 'BC_split' / f'{sample}.fastq.bz2'),
                'umi_trimmed_fastq_MB': self._file_size_mb(
                    data_dir / 'UMI_trimmed' / f'{sample}_UMI-trimmed.fastq.bz2'),
                'alignment_json_MB': self._file_size_mb(
                    data_dir / 'SWalign' / f'{sample}_SWalign.json.bz2'),
                'per_read_stats_MB': self._file_size_mb(
                    data_dir / 'stats_collection' / f'{sample}_stats.csv.bz2'),
            }

    def _save_metrics(self):
        """Save computing_metrics.csv (or per-sample JSON in single-sample mode)."""
        if self.sample_index is not None:
            # Single-sample SLURM mode: save per-sample JSON
            metrics_dir = self.project_dir / 'data' / 'metrics'
            metrics_dir.mkdir(exist_ok=True)

            sample_name = self.sample_df.iloc[self.sample_index]['sample_name_unique']
            sample_metrics = {
                'sample_name_unique': sample_name,
                **self.file_metrics.get(sample_name, {}),
                **{f'{k}_sec': round(v, 2) for k, v in self.stage_timings.items()},
            }

            out_path = metrics_dir / f'{sample_name}_metrics.json'
            with open(out_path, 'w') as f:
                json.dump(sample_metrics, f, indent=2)
            self.log(f"  Saved per-sample metrics: {out_path}")
        else:
            # Full pipeline or aggregation mode: save CSV
            self._aggregate_metrics_csv()

    def _aggregate_metrics_csv(self):
        """Aggregate per-sample JSONs (if any) and pipeline metrics into CSV."""
        metrics_dir = self.project_dir / 'data' / 'metrics'
        rows = []

        # Try loading per-sample JSONs from SLURM array jobs
        if metrics_dir.exists():
            for json_file in sorted(metrics_dir.glob('*_metrics.json')):
                with open(json_file) as f:
                    rows.append(json.load(f))

        if rows:
            # Merge with any file metrics collected in this run
            df = pd.DataFrame(rows)
        else:
            # No per-sample JSONs — build from current run's data
            self._collect_file_metrics()
            samples = self.sample_df['sample_name_unique'].tolist()
            for sample in samples:
                row = {'sample_name_unique': sample}
                row.update(self.file_metrics.get(sample, {}))
                rows.append(row)
            df = pd.DataFrame(rows)

        # Add pipeline-level timings as columns (same value for all rows)
        for stage, duration in self.stage_timings.items():
            col = f'{stage}_sec'
            if col not in df.columns:
                df[col] = round(duration, 2)

        # Add aggregate stats file size
        agg_stats = self.project_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
        if agg_stats.exists():
            df['aggregate_stats_MB'] = self._file_size_mb(agg_stats)

        out_path = self._ensure_dir('logs') / 'computing_metrics.csv'
        df.to_csv(out_path, index=False)
        self.log(f"  Saved computing metrics: {out_path} ({len(df)} samples)")

    def save_outputs(self):
        """Save output files"""
        self.log("=" * 60)
        self.log("Saving outputs")
        self.log("=" * 60)

        # Save inp_file_df
        inp_file_output = self.project_dir / 'inp_file_df.xlsx'
        self.inp_file_df.to_excel(inp_file_output, index=False)
        self.log(f"  Saved: {inp_file_output}")

        # Save sample_df
        sample_df_output = self.project_dir / 'sample_df.xlsx'
        self.sample_df.to_excel(sample_df_output, index=False)
        self.log(f"  Saved: {sample_df_output}")

        # ALL_stats_aggregate.csv lives in data/stats_collection/ (no root copy)
        stats_file = self.project_dir / 'data' / 'stats_collection' / 'ALL_stats_aggregate.csv'
        if stats_file.exists():
            self.log(f"  Stats: data/stats_collection/ALL_stats_aggregate.csv")
        else:
            self.log(f"  WARNING: ALL_stats_aggregate.csv not found", level="WARN")

        # Charge results live in results/charge/ (no root copies)
        if hasattr(self, 'charge_results') and self.charge_results:
            for level in self.charge_results:
                self.log(f"  Charge: results/charge/charge_df_{level}.csv")

    @staticmethod
    def parse_stages(stages_str):
        """
        Parse --stages argument into a set of stage identifiers.

        Accepts:
            'all'       -> all stages
            '0-2'       -> stages 0a, 0b, 0c, 1, 2
            '3-5'       -> stages 3, 4, 5
            '0a,1,3'    -> specific stages
            '2,3,5'     -> specific stages

        Returns:
            Set of stage identifiers: {'0a', '0b', '0c', '1', '2', '3', '4', '5'}
        """
        ALL_STAGES = {'0a', '0b', '0c', '1', '2', '3', '4', '5', '6', '7'}

        if stages_str == 'all':
            return ALL_STAGES

        stages = set()
        for part in stages_str.split(','):
            part = part.strip()
            if '-' in part and not part.startswith('0'):
                # Range like '3-5'
                start, end = part.split('-')
                # Expand numeric ranges
                for i in range(int(start), int(end) + 1):
                    stages.add(str(i))
            elif part == '0-2' or part == '0-1':
                # Special: includes sub-stages
                stages.update({'0a', '0b', '0c'})
                if part == '0-2':
                    stages.update({'1', '2'})
                else:
                    stages.add('1')
            elif part in ALL_STAGES:
                stages.add(part)
            elif part == '0':
                stages.update({'0a', '0b', '0c'})
            else:
                raise ValueError(f"Unknown stage: '{part}'. Valid: {sorted(ALL_STAGES)}")

        return stages

    def _load_prerequisites(self, stages):
        """
        Load data prerequisites when skipping early stages.

        When stages 0-1 are skipped, we still need sample_df and tRNA_data
        for later stages. Load them from existing files.
        """
        # Always need sample_df -- load from existing file or from config
        skip_early = not stages.intersection({'0a', '0b', '0c', '1'})
        if skip_early:
            # Try loading existing sample_df from output directory
            # Check both output root and data/ subdirectory
            sample_df_path = self.project_dir / 'sample_df.xlsx'
            if not sample_df_path.exists():
                sample_df_path = self.project_dir / 'data' / 'sample_df.xlsx'

            if sample_df_path.exists():
                self.log(f"Loading existing sample_df from {sample_df_path}")
                self.sample_df = pd.read_excel(sample_df_path)
                # Also load inp_file_df if it exists
                inp_file_path = self.project_dir / 'inp_file_df.xlsx'
                if not inp_file_path.exists():
                    inp_file_path = self.project_dir / 'data' / 'inp_file_df.xlsx'
                if inp_file_path.exists():
                    self.inp_file_df = pd.read_excel(inp_file_path)
                else:
                    self.inp_file_df = pd.DataFrame()
            else:
                # Fall back to loading from config
                self.log("No existing sample_df found, loading from config...")
                self.load_sample_info()

        # Stage 2 needs tRNA_data -- load from config if stage 1 was skipped
        if '2' in stages and '1' not in stages:
            if not hasattr(self, 'tRNA_data') or self.tRNA_data is None:
                self.log("Loading tRNA database for stats collection...")
                tRNA_database = self.config['tRNA_database']
                self.tRNA_data = read_tRNAdb_info(tRNA_database)

    def run(self, stages=None):
        """
        Run preprocessing pipeline.

        Args:
            stages: Set of stage identifiers to run, or None for all stages.
                    e.g. {'2', '3', '5'} to run only stats, charge, and QC.
        """
        if stages is None:
            stages = {'0a', '0b', '0c', '1', '2', '3', '4', '5', '6', '7'}

        try:
            self.log("Starting tRNA-charge-seq preprocessing pipeline")
            self.log(f"Configuration: {self.config_file}")
            self.log(f"Project directory: {self.project_dir}")
            self.log(f"Parallel jobs: {self.n_jobs}")
            self.log(f"CPU budget: {self.n_jobs} jobs × {self.threads_per_job} threads = "
                     f"{self.n_jobs * self.threads_per_job} cores")
            self.log(f"Stages to run: {', '.join(sorted(stages))}")

            if self.sample_index is not None:
                self.log(f"Sample index: {self.sample_index} (single-sample mode)")
                if stages.intersection({'0a', '0b'}):
                    self.log("WARNING: --sample-index has no effect on stages 0a/0b "
                             "(they process all file pairs)", level="WARN")

            # Validate project directory for early stages
            if stages.intersection({'0a', '0b', '0c', '1'}):
                seq_dir = self.config.get('seq_dir', 'raw_fastq')
                raw_data_path = self.project_dir / 'data' / seq_dir
                if not raw_data_path.is_dir():
                    raise FileNotFoundError(
                        f"Raw data directory not found: {raw_data_path}\n"
                        f"--project-dir must point to the project root containing data/{seq_dir}/"
                    )

            # Setup directories
            self.setup_directories()

            # Load sample info for early stages
            if stages.intersection({'0a', '0b'}):
                # Stages 0a/0b: load fresh from config (sample list Excel)
                self.load_sample_info()
            elif stages.intersection({'0c', '1'}):
                # Stages 0c/1 without 0a/0b: need sample_df with N_total from 0b
                sample_df_path = self.project_dir / 'sample_df.xlsx'
                if sample_df_path.exists():
                    self.log(f"Loading existing sample_df from {sample_df_path}")
                    self.sample_df = pd.read_excel(sample_df_path)
                    inp_file_path = self.project_dir / 'inp_file_df.xlsx'
                    if inp_file_path.exists():
                        self.inp_file_df = pd.read_excel(inp_file_path)
                    else:
                        self.inp_file_df = pd.DataFrame()
                    self.log(f"Loaded {len(self.sample_df)} samples")
                    self.log(f"Input files: {len(self.inp_file_df)}")
                else:
                    self.log("WARNING: No sample_df.xlsx found — loading from config. "
                             "Stage 0c may fail if N_total column is missing.", level="WARN")
                    self.load_sample_info()

            # Load prerequisites for later stages when skipping early ones
            self._load_prerequisites(stages)

            # Run preprocessing stages (0-2)
            if '0a' in stages:
                t0 = time.time()
                self.stage_0a_adapter_removal()
                self.stage_timings['stage_0a'] = time.time() - t0
            if '0b' in stages:
                t0 = time.time()
                self.stage_0b_barcode_split()
                self.stage_timings['stage_0b'] = time.time() - t0
            if '0c' in stages:
                t0 = time.time()
                self.stage_0c_umi_trim()
                self.stage_timings['stage_0c'] = time.time() - t0
            if '1' in stages:
                t0 = time.time()
                self.stage_1_alignment()
                self.stage_timings['stage_1'] = time.time() - t0
            if '2' in stages:
                t0 = time.time()
                self.stage_2_stats_collection()
                self.stage_timings['stage_2'] = time.time() - t0

            # Run analysis stages (3-5)
            if '3' in stages:
                if self.config.get('run_charge_quantification', True):
                    t0 = time.time()
                    self.stage_3_charge_quantification()
                    self.stage_timings['stage_3'] = time.time() - t0
                else:
                    self.log("Skipping Stage 3: Charge Quantification (disabled in config)")

                # Fragment analysis (part of stage 3)
                if self.config.get('run_fragment_analysis', True):
                    t0 = time.time()
                    self._run_fragment_analysis()
                    self.stage_timings['stage_3_fragments'] = time.time() - t0
                else:
                    self.log("Skipping fragment analysis (disabled in config)")

            if '4' in stages:
                if self.config.get('run_parquet_storage', False):
                    t0 = time.time()
                    self.stage_4_parquet_storage()
                    self.stage_timings['stage_4'] = time.time() - t0
                else:
                    self.log("Skipping Stage 4: Parquet Storage (disabled in config)")

            if '5' in stages:
                if self.config.get('run_qc_report', True):
                    t0 = time.time()
                    self.stage_5_qc_report()
                    self.stage_timings['stage_5'] = time.time() - t0
                else:
                    self.log("Skipping Stage 5: QC Report (disabled in config)")

            if '6' in stages:
                if self.config.get('run_modification_analysis', False):
                    t0 = time.time()
                    self.stage_6_modification_analysis()
                    self.stage_timings['stage_6'] = time.time() - t0
                else:
                    self.log("Skipping Stage 6: Modification Analysis (disabled in config)")

            if '7' in stages:
                if self.config.get('run_abundance_analysis', False):
                    t0 = time.time()
                    self.stage_7_abundance_analysis()
                    self.stage_timings['stage_7'] = time.time() - t0
                else:
                    self.log("Skipping Stage 7: Abundance Analysis (disabled in config)")

            # Collect file metrics and save
            if self.sample_index is not None:
                self._collect_file_metrics()
                self._save_metrics()
            else:
                self._collect_file_metrics()

            # Save outputs (skip in single-sample mode)
            if self.sample_index is None:
                self.save_outputs()
                self._save_metrics()
            else:
                self.log("Skipping save_outputs in single-sample mode")

            # Final summary
            self.status['end_time'] = datetime.now()
            duration = self.status['end_time'] - self.status['start_time']

            self.log("=" * 60)
            self.log("Pipeline Complete!")
            self.log("=" * 60)
            self.log(f"Duration: {duration}")
            self.log(f"Stages completed: {', '.join(self.status['stages_completed'])}")
            if self.stage_timings:
                self.log(f"Stage timings:")
                for stage, secs in sorted(self.stage_timings.items()):
                    m, s = divmod(int(secs), 60)
                    self.log(f"  {stage}: {m}m {s}s")
            self.log(f"Output files:")
            self.log(f"  - inp_file_df.xlsx")
            self.log(f"  - sample_df.xlsx")
            self.log(f"  - data/stats_collection/ALL_stats_aggregate.csv")

            # List charge results if generated
            if '3' in self.status['stages_completed']:
                self.log(f"  - results/charge/charge_df_*.csv")
                self.log(f"  - results/charge/charge_summary.csv")
                if hasattr(self, '_fragment_analysis_ran'):
                    self.log(f"  - results/fragments/*.parquet")

            if '5' in self.status['stages_completed']:
                self.log(f"  - qc_reports/QC_summary.csv")
                self.log(f"  - qc_reports/QC_report.html")

            self.log(f"  - logs/computing_metrics.csv")
            self.log("")
            if '3' not in self.status['stages_completed']:
                self.log("Next step: Run charge quantification (Stage 3)")
                self.log("  Re-run with --stages 3 or:")
                self.log("  python -m trnaseq quantify -i ALL_stats_aggregate.csv -o charge.csv")

        except Exception as e:
            self.log(f"Pipeline failed: {str(e)}", level="ERROR")
            raise


def main():
    parser = argparse.ArgumentParser(
        description="Unified tRNA-charge-seq preprocessing pipeline (Stages 0-7)"
    )
    parser.add_argument(
        '--config', required=True,
        help='YAML configuration file'
    )
    parser.add_argument(
        '--project-dir', required=False, default=None,
        help='Project directory (must contain data/raw_fastq/ for stages 0-1)'
    )
    parser.add_argument(
        '--output-dir', required=False, default=None, dest='output_dir_deprecated',
        help='Deprecated: use --project-dir instead'
    )
    parser.add_argument(
        '--n-jobs', type=int, default=4,
        help='Number of parallel jobs (default: 4)'
    )
    parser.add_argument(
        '--skip-charge', action='store_true',
        help='Skip charge quantification (Stage 3)'
    )
    parser.add_argument(
        '--parquet', action='store_true',
        help='Save to Parquet format (Stage 4)'
    )
    parser.add_argument(
        '--stages', type=str, default='all',
        help='Stages to run: "all", "0-2", "3-5", or specific like "0a,0b,1,3"'
    )
    parser.add_argument(
        '--sample-index', type=int, default=None,
        help='Process only sample at this 0-based index (for SLURM array jobs)'
    )
    parser.add_argument(
        '--threads-per-job', type=int, default=None,
        help='Threads per subprocess (AdapterRemoval/SWIPE). '
             'Default: config threads_per_job or 2'
    )
    parser.add_argument(
        '--preflight', action='store_true',
        help='Validate config, files, database, and tools without running '
             'the pipeline. Exits with code 0 (pass) or 1 (fail).'
    )

    args = parser.parse_args()

    # Resolve --project-dir / --output-dir (deprecated)
    project_dir = args.project_dir or args.output_dir_deprecated
    if project_dir is None:
        parser.error('--project-dir is required')
    if args.output_dir_deprecated and not args.project_dir:
        print('Warning: --output-dir is deprecated, use --project-dir instead',
              file=sys.stderr)

    # Run pipeline
    pipeline = PreprocessingPipeline(
        config_file=args.config,
        project_dir=project_dir,
        n_jobs=args.n_jobs,
        sample_index=args.sample_index,
        threads_per_job=args.threads_per_job,
    )

    # Override config with CLI arguments
    if args.skip_charge:
        pipeline.config['run_charge_quantification'] = False
    if args.parquet:
        pipeline.config['run_parquet_storage'] = True

    # Parse stages
    stages = PreprocessingPipeline.parse_stages(args.stages)

    # Preflight mode: validate and exit
    if args.preflight:
        sys.exit(pipeline.preflight(stages=stages))

    pipeline.run(stages=stages)


if __name__ == '__main__':
    main()
