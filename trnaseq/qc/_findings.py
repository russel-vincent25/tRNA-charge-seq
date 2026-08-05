"""
Findings registry — the threshold-to-sentence layer
====================================================

Each panel computes a metric. A :class:`FindingRule` knows how to turn a
metric value into a one-sentence diagnostic. :func:`evaluate_findings`
runs every rule registered for a panel anchor and yields :class:`Finding`
records the panel can render under itself (T2) and the index page can
aggregate into a run-health roll-up (T1).

Adding a rule
-------------
``register(FindingRule(...))`` appends to the global :data:`FINDINGS`. The
seed rules at the bottom of this module are the ones called out in §4.3 of
the report-improvements PRP; later tasks extend the list.

Severity has three levels: ``pass`` (no problem worth flagging), ``warn``
(investigate), ``fail`` (sample should be quarantined). Each panel renders
the **highest severity** finding it triggered, and the index page rolls up
the worst severity across the whole run.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable, Mapping, Optional

import pandas as pd


# ---------------------------------------------------------------------------
# Severity ordering — used everywhere we sort findings.
# ---------------------------------------------------------------------------

#: Increasing severity. Index into this list yields the comparison rank.
SEVERITY_ORDER: tuple[str, ...] = ('pass', 'warn', 'fail')


def severity_rank(severity: str) -> int:
    """Numeric rank: pass=0, warn=1, fail=2. Unknown → -1."""
    try:
        return SEVERITY_ORDER.index(severity)
    except ValueError:
        return -1


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FindingRule:
    """Predicate + template for one diagnostic.

    Parameters
    ----------
    anchor : str
        Panel anchor this rule applies to. Conventionally ``"<report>-<panel>"``,
        matching the ``anchor=`` kwarg passed to
        :func:`trnaseq.qc._common.render_panel`. Wildcards aren't supported
        — register a rule per anchor (or share a predicate).
    predicate : Callable[[float], bool]
        Returns True when the rule should fire. Receives the metric value.
    severity : str
        One of ``pass``, ``warn``, ``fail``. Higher overrides lower.
    template : str
        ``str.format``-compatible template. The full row dict is passed in
        as keyword args, so ``{sample}``, ``{value}``, ``{aa}`` etc. all
        work — whatever columns the metric DataFrame carries.
    threshold : Any, optional
        Bookkeeping field; surfaces in the template as ``{threshold}``.
    """
    anchor: str
    predicate: Callable[[float], bool]
    severity: str
    template: str
    threshold: Any = None


@dataclass(frozen=True)
class Finding:
    """A triggered rule — what gets rendered/aggregated."""
    anchor: str
    severity: str
    message: str
    sample: Optional[str] = None
    value: Optional[float] = None


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

#: Module-level registry. Mutated by :func:`register`.
FINDINGS: list[FindingRule] = []


def register(rule: FindingRule) -> FindingRule:
    """Add a rule to :data:`FINDINGS` and return it (handy as a decorator-ish call)."""
    FINDINGS.append(rule)
    return rule


def rules_for(anchor: str) -> list[FindingRule]:
    """Return every rule registered for ``anchor`` (order-preserving)."""
    return [r for r in FINDINGS if r.anchor == anchor]


def clear() -> None:
    """Drop all registered rules. Used by tests to start from a clean slate."""
    FINDINGS.clear()


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate_findings(
    anchor: str,
    metric_df: pd.DataFrame,
    *,
    value_col: str = 'value',
) -> list[Finding]:
    """Apply all rules for ``anchor`` to every row of ``metric_df``.

    Each row may trigger 0–N rules; if multiple rules fire on the same row
    we keep only the **highest severity** one (since we render at most one
    sentence per row anyway). Across rows, returns the full list.

    Parameters
    ----------
    anchor : str
        Panel anchor the rules are registered under.
    metric_df : pd.DataFrame
        Must contain ``value_col``; any other columns become template
        substitutions. NaN values in ``value_col`` are silently skipped.
    value_col : str
        Name of the metric column. Defaults to ``'value'``.
    """
    if not isinstance(metric_df, pd.DataFrame):
        raise TypeError(f'metric_df must be a DataFrame, got {type(metric_df).__name__}')
    rules = rules_for(anchor)
    if not rules or metric_df.empty or value_col not in metric_df.columns:
        return []

    out: list[Finding] = []
    for _, row in metric_df.iterrows():
        v = row[value_col]
        if v is None or (isinstance(v, float) and pd.isna(v)):
            continue
        row_dict = row.to_dict()
        # Make {value} always available even if row has no 'value' column.
        row_dict.setdefault('value', v)
        triggered: list[tuple[int, Finding]] = []
        for rule in rules:
            try:
                fires = bool(rule.predicate(v))
            except (TypeError, ValueError):
                fires = False
            if not fires:
                continue
            fmt = dict(row_dict)
            fmt.setdefault('threshold', rule.threshold)
            try:
                msg = rule.template.format(**fmt)
            except (KeyError, IndexError, ValueError):
                msg = rule.template
            triggered.append((
                severity_rank(rule.severity),
                Finding(
                    anchor=anchor,
                    severity=rule.severity,
                    message=msg,
                    sample=str(row_dict.get('sample')) if 'sample' in row_dict else None,
                    value=float(v) if isinstance(v, (int, float)) else None,
                ),
            ))
        if triggered:
            triggered.sort(key=lambda t: t[0], reverse=True)
            out.append(triggered[0][1])
    return out


def evaluate_value(
    anchor: str,
    value: float,
    *,
    sample: Optional[str] = None,
    extra: Optional[Mapping[str, Any]] = None,
) -> Optional[Finding]:
    """Single-value convenience: highest-severity finding for one metric.

    Returns ``None`` if no rule fires.
    """
    rules = rules_for(anchor)
    if not rules or value is None:
        return None
    row: dict[str, Any] = {'value': value, 'sample': sample}
    if extra:
        row.update(extra)
    triggered: list[tuple[int, Finding]] = []
    for rule in rules:
        try:
            fires = bool(rule.predicate(value))
        except (TypeError, ValueError):
            fires = False
        if not fires:
            continue
        fmt = dict(row)
        fmt.setdefault('threshold', rule.threshold)
        try:
            msg = rule.template.format(**fmt)
        except (KeyError, IndexError, ValueError):
            msg = rule.template
        triggered.append((
            severity_rank(rule.severity),
            Finding(
                anchor=anchor,
                severity=rule.severity,
                message=msg,
                sample=sample,
                value=float(value) if isinstance(value, (int, float)) else None,
            ),
        ))
    if not triggered:
        return None
    triggered.sort(key=lambda t: t[0], reverse=True)
    return triggered[0][1]


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate_findings_by_severity(findings: Iterable[Finding]) -> dict[str, list[Finding]]:
    """Group findings by ``severity`` — keys: 'pass', 'warn', 'fail'."""
    out: dict[str, list[Finding]] = {sev: [] for sev in SEVERITY_ORDER}
    for f in findings:
        out.setdefault(f.severity, []).append(f)
    return out


def worst_severity(findings: Iterable[Finding]) -> str:
    """Return the highest severity present, or ``'pass'`` if empty/clean."""
    worst = 'pass'
    worst_rank = severity_rank(worst)
    for f in findings:
        r = severity_rank(f.severity)
        if r > worst_rank:
            worst, worst_rank = f.severity, r
    return worst


def summarize_findings(findings: Iterable[Finding]) -> dict[str, Any]:
    """Counts + worst severity, ready for the index-page roll-up (T1)."""
    grouped = aggregate_findings_by_severity(findings)
    return {
        'n_pass': len(grouped.get('pass', [])),
        'n_warn': len(grouped.get('warn', [])),
        'n_fail': len(grouped.get('fail', [])),
        'worst': worst_severity(findings),
    }


# ---------------------------------------------------------------------------
# Seed rules (PRP §4.3) — extended by later tasks.
# ---------------------------------------------------------------------------

# QC report — SWIPE mapping rate per sample
register(FindingRule(
    anchor='qc-mapping',
    predicate=lambda x: x < 40,
    severity='fail',
    template='Sample {sample} mapping {value:.1f}% (fail threshold 40%)',
    threshold=40,
))
register(FindingRule(
    anchor='qc-mapping',
    predicate=lambda x: 40 <= x < 70,
    severity='warn',
    template='Sample {sample} mapping {value:.1f}% (warn threshold 70%)',
    threshold=70,
))

# Charge report — mean charging fraction per sample
register(FindingRule(
    anchor='charge-mean',
    predicate=lambda x: x < 30,
    severity='fail',
    template='Sample {sample} mean charge {value:.1f}% — investigate RNA quality',
    threshold=30,
))
register(FindingRule(
    anchor='charge-mean',
    predicate=lambda x: 30 <= x < 50,
    severity='warn',
    template='Sample {sample} mean charge {value:.1f}% (warn threshold 50%)',
    threshold=50,
))

# Charge report — per-amino-acid charging
register(FindingRule(
    anchor='charge-aa',
    predicate=lambda x: x < 20,
    severity='warn',
    template='Charge for {aa} in {sample} is {value:.1f}% — possible synthetase issue',
    threshold=20,
))


__all__ = [
    'FindingRule',
    'Finding',
    'FINDINGS',
    'SEVERITY_ORDER',
    'register',
    'rules_for',
    'clear',
    'severity_rank',
    'evaluate_findings',
    'evaluate_value',
    'aggregate_findings_by_severity',
    'worst_severity',
    'summarize_findings',
]
