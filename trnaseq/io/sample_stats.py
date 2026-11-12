"""
Per-sample stage 0c/1 stats files.

SLURM array tasks (``--sample-index``) run stages 0c and 1 for one sample and
exit without writing ``sample_df.xlsx``; hundreds of tasks would otherwise race
on one spreadsheet. Each task instead writes that stage's columns for its
sample to ``data/sample_stats/{sample}.stage{stage}.json``. The aggregation run
(stages 2-7) merges these files into ``sample_df`` before stage 2, so the QC
report, the stage 2 ``Mapping_percent`` and ``sample_df.xlsx`` get the same
values a single-process run would produce.

Files are written atomically (temp file + ``os.replace``), so a reader never
sees a half-written file, and their values overwrite whatever ``sample_df``
already holds for that sample, because the file is always the newer record.
"""

from __future__ import annotations

import json
import math
import os
import socket
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SCHEMA_VERSION = 1
SAMPLE_STATS_SUBDIR = Path('data') / 'sample_stats'

#: Columns each stage adds to sample_df.
STAGE_COLUMNS: Dict[str, Tuple[str, ...]] = {
    '0c': ('N_after_trim', 'N_UMI_observed', 'N_UMI_expected', 'N_after_downsample',
           'percent_seqs_after_UMI_trim', 'percent_UMI_obs-vs-exp'),
    '1': ('N_mapped', 'percent_single_annotation', 'percent_multiple_annotation',
          'percent_multiple_codons', 'Mapping_percent'),
}

#: A sample already has a stage's stats when this column is non-null in
#: sample_df (e.g. a sample_df.xlsx written by a single-process run). Only the
#: stage itself writes these columns; stage 2 rewrites Mapping_percent, so that
#: one cannot serve as the marker for stage 1.
COVERAGE_COLUMN: Dict[str, str] = {'0c': 'N_after_trim', '1': 'N_mapped'}


@dataclass
class StageStatus:
    """Outcome of merging one stage's stats files."""
    stage: str
    n_samples: int                      # samples in scope
    n_files: int                        # stats files found and merged
    missing: List[Tuple[int, str]] = field(default_factory=list)
    # (row position in sample_df == SLURM array index, sample_name_unique)
    # for samples with neither a file nor values already in sample_df

    @property
    def complete(self) -> bool:
        return not self.missing


# ---------------------------------------------------------------------------
# Paths and I/O
# ---------------------------------------------------------------------------

def sample_stats_dir(project_dir) -> Path:
    return Path(project_dir) / SAMPLE_STATS_SUBDIR


def sample_stats_path(project_dir, sample: str, stage: str) -> Path:
    return sample_stats_dir(project_dir) / f'{sample}.stage{stage}.json'


def _to_json(value):
    """Convert a pandas/numpy scalar to a JSON-safe value (NaN -> null)."""
    if isinstance(value, np.generic):
        value = value.item()
    if value is None or value is pd.NA:
        return None
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _provenance() -> dict:
    env = os.environ
    return {
        'written_at': datetime.now().isoformat(timespec='seconds'),
        'host': socket.gethostname(),
        'slurm_job_id': env.get('SLURM_JOB_ID'),
        'slurm_array_job_id': env.get('SLURM_ARRAY_JOB_ID'),
        'slurm_array_task_id': env.get('SLURM_ARRAY_TASK_ID'),
    }


def _atomic_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    try:
        with open(tmp, 'w') as fh:
            json.dump(payload, fh, indent=2)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def write_sample_stats(project_dir, result_df: pd.DataFrame, stage: str) -> List[Path]:
    """Write one stats file per row of ``result_df`` for ``stage``.

    Only the columns in ``STAGE_COLUMNS[stage]`` that ``result_df`` has are
    written. Returns the paths written.
    """
    if stage not in STAGE_COLUMNS:
        raise ValueError(f"Unknown stage {stage!r}; expected one of {sorted(STAGE_COLUMNS)}")
    cols = [c for c in STAGE_COLUMNS[stage] if c in result_df.columns]
    if not cols:
        return []
    prov = _provenance()
    paths = []
    for _, row in result_df.iterrows():
        sample = row['sample_name_unique']
        payload = {
            'schema': SCHEMA_VERSION,
            'sample_name_unique': sample,
            'stage': stage,
            'values': {c: _to_json(row[c]) for c in cols},
            'provenance': prov,
        }
        path = sample_stats_path(project_dir, sample, stage)
        _atomic_write_json(path, payload)
        paths.append(path)
    return paths


def read_sample_stats(path: Path, sample: str, stage: str) -> dict:
    """Read one stats file and check it belongs to ``sample`` and ``stage``."""
    try:
        with open(path) as fh:
            payload = json.load(fh)
    except json.JSONDecodeError as err:
        raise ValueError(f"{path}: not valid JSON ({err})") from err
    if payload.get('schema', 0) > SCHEMA_VERSION:
        raise ValueError(f"{path}: schema {payload.get('schema')} is newer than this "
                         f"trnaseq supports ({SCHEMA_VERSION}); update the repo")
    if payload.get('sample_name_unique') != sample or payload.get('stage') != stage:
        raise ValueError(f"{path}: holds sample {payload.get('sample_name_unique')!r} "
                         f"stage {payload.get('stage')!r}, expected {sample!r} stage {stage!r}")
    return payload


# ---------------------------------------------------------------------------
# Merging into sample_df
# ---------------------------------------------------------------------------

def assign_sample_values(sample_df: pd.DataFrame, values: pd.DataFrame) -> pd.DataFrame:
    """Overwrite ``sample_df`` cells with ``values`` (indexed by sample_name_unique).

    Adds columns as needed. Row order and every cell not named in ``values``
    are left unchanged; samples absent from ``values`` keep their old values.
    """
    if values.empty:
        return sample_df
    if not values.index.is_unique:
        dup = values.index[values.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate samples in values: {dup}")
    df = sample_df.copy()
    has = df['sample_name_unique'].isin(values.index)
    for col in values.columns:
        mapped = df['sample_name_unique'].map(values[col])
        df[col] = mapped.where(has, df[col]) if col in df.columns else mapped
    return df


def merge_sample_stats(sample_df: pd.DataFrame, project_dir,
                       stages: Sequence[str] = ('0c', '1'),
                       samples: Optional[Iterable[str]] = None,
                       ) -> Tuple[pd.DataFrame, Dict[str, StageStatus]]:
    """Merge per-sample stats files into ``sample_df``.

    Parameters:
        sample_df: Sample table; its row order defines the SLURM array indices.
        project_dir: Project root holding ``data/sample_stats/``.
        stages: Stages whose files to merge.
        samples: Restrict to these samples (default: every sample in sample_df).

    Returns:
        (merged sample_df, {stage: StageStatus})
    """
    df = sample_df
    names = df['sample_name_unique']
    scope = set(names) if samples is None else set(samples)
    in_scope = names.isin(scope)
    status = {}
    for stage in stages:
        records = {}
        for sample in names[in_scope]:
            path = sample_stats_path(project_dir, sample, stage)
            if path.exists():
                records[sample] = read_sample_stats(path, sample, stage)['values']
        if records:
            values = (pd.DataFrame.from_dict(records, orient='index')
                      .apply(pd.to_numeric, errors='coerce'))
            df = assign_sample_values(df, values)

        cov = COVERAGE_COLUMN[stage]
        missing_mask = in_scope & (df[cov].isna() if cov in df.columns else True)
        positions = np.flatnonzero(np.asarray(missing_mask, dtype=bool))
        status[stage] = StageStatus(
            stage=stage,
            n_samples=int(in_scope.sum()),
            n_files=len(records),
            missing=[(int(i), str(df['sample_name_unique'].iat[i])) for i in positions],
        )
    return df, status


def format_array_indices(indices: Iterable[int]) -> str:
    """Compress indices to sbatch ``--array`` syntax: [0, 1, 2, 5] -> '0-2,5'."""
    idx = sorted(set(int(i) for i in indices))
    parts = []
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
            j += 1
        parts.append(str(idx[i]) if i == j else f'{idx[i]}-{idx[j]}')
        i = j + 1
    return ','.join(parts)


# ---------------------------------------------------------------------------
# Derived columns
# ---------------------------------------------------------------------------

def recompute_mapping_percent(sample_df: pd.DataFrame) -> pd.DataFrame:
    """Set ``Mapping_percent = N_total_aligned / reads aligned * 100``.

    The denominator is ``N_after_downsample``, the number of reads SWIPE was
    given; ``N_after_trim`` is only used for rows that lack it. Dividing by
    ``N_after_trim`` when downsampling is on would understate the mapping rate
    by N_after_downsample / N_after_trim. A zero denominator gives 0. Rows with
    neither column keep their existing ``Mapping_percent`` (from stage 1).
    """
    if 'N_total_aligned' not in sample_df.columns:
        return sample_df
    denom = pd.Series(np.nan, index=sample_df.index)
    for col in ('N_after_downsample', 'N_after_trim'):
        if col in sample_df.columns:
            denom = denom.fillna(pd.to_numeric(sample_df[col], errors='coerce'))
    known = denom.notna()
    if not known.any():
        return sample_df
    aligned = pd.to_numeric(sample_df['N_total_aligned'], errors='coerce')
    new = (aligned / denom.where(denom > 0) * 100).where(denom > 0, 0.0)
    df = sample_df.copy()
    if 'Mapping_percent' in df.columns:
        df['Mapping_percent'] = new.where(known, df['Mapping_percent'])
    else:
        df['Mapping_percent'] = new.where(known)
    return df
