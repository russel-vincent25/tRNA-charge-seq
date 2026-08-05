"""
Tests for trnaseq.qc._findings — threshold → diagnostic sentence layer.
"""

from __future__ import annotations

import pandas as pd
import pytest

from trnaseq.qc import _findings
from trnaseq.qc._findings import (
    FINDINGS,
    FindingRule,
    SEVERITY_ORDER,
    aggregate_findings_by_severity,
    evaluate_findings,
    evaluate_value,
    register,
    rules_for,
    severity_rank,
    summarize_findings,
    worst_severity,
)


@pytest.fixture
def isolated_registry():
    """Snapshot FINDINGS, swap it for an empty list for the test, restore on exit."""
    saved = list(FINDINGS)
    FINDINGS.clear()
    try:
        yield FINDINGS
    finally:
        FINDINGS.clear()
        FINDINGS.extend(saved)


# ---------------------------------------------------------------------------
# Severity utilities
# ---------------------------------------------------------------------------

class TestSeverity:
    def test_pass_warn_fail_ordering(self):
        assert SEVERITY_ORDER == ('pass', 'warn', 'fail')

    def test_severity_rank_monotonic(self):
        assert severity_rank('pass') < severity_rank('warn') < severity_rank('fail')

    def test_unknown_severity_returns_negative(self):
        assert severity_rank('catastrophic') == -1

    def test_worst_severity_picks_max(self):
        findings = [
            _findings.Finding(anchor='x', severity='warn', message='m'),
            _findings.Finding(anchor='x', severity='fail', message='m'),
            _findings.Finding(anchor='x', severity='pass', message='m'),
        ]
        assert worst_severity(findings) == 'fail'

    def test_worst_severity_clean_run_is_pass(self):
        assert worst_severity([]) == 'pass'


# ---------------------------------------------------------------------------
# Seed rules — PRP §4.3
# ---------------------------------------------------------------------------

class TestSeedRules:
    """The seed rules from PRP §4.3 must be registered on import."""

    def test_qc_mapping_fail_rule_present(self):
        rules = rules_for('qc-mapping')
        fails = [r for r in rules if r.severity == 'fail']
        assert len(fails) >= 1
        assert any(r.predicate(30) for r in fails)         # fires at 30%
        assert not any(r.predicate(40) for r in fails)     # not at 40

    def test_qc_mapping_warn_rule_present(self):
        rules = rules_for('qc-mapping')
        warns = [r for r in rules if r.severity == 'warn']
        assert len(warns) >= 1
        assert any(r.predicate(60) for r in warns)         # 60 is in warn band
        assert not any(r.predicate(80) for r in warns)     # 80 is clean

    def test_charge_mean_fail_rule_present(self):
        rules = rules_for('charge-mean')
        fails = [r for r in rules if r.severity == 'fail']
        assert any(r.predicate(20) for r in fails)
        assert not any(r.predicate(50) for r in fails)

    def test_charge_aa_warn_rule_present(self):
        rules = rules_for('charge-aa')
        warns = [r for r in rules if r.severity == 'warn']
        assert any(r.predicate(10) for r in warns)


# ---------------------------------------------------------------------------
# Registry operations (use isolated_registry to keep seeds out of the way)
# ---------------------------------------------------------------------------

class TestRegistry:
    def test_register_appends(self, isolated_registry):
        r = FindingRule('x', lambda v: True, 'warn', 'always')
        register(r)
        assert FINDINGS == [r]

    def test_rules_for_filters(self, isolated_registry):
        a = FindingRule('a', lambda v: True, 'warn', 'A')
        b = FindingRule('b', lambda v: True, 'warn', 'B')
        c = FindingRule('a', lambda v: True, 'fail', 'A2')
        for r in (a, b, c):
            register(r)
        out = rules_for('a')
        assert out == [a, c]

    def test_clear_empties_registry(self, isolated_registry):
        register(FindingRule('x', lambda v: True, 'warn', 'tpl'))
        assert FINDINGS
        _findings.clear()
        assert FINDINGS == []


# ---------------------------------------------------------------------------
# evaluate_findings
# ---------------------------------------------------------------------------

class TestEvaluateFindings:
    def test_no_rules_returns_empty(self, isolated_registry):
        df = pd.DataFrame({'value': [10, 20], 'sample': ['a', 'b']})
        assert evaluate_findings('whatever', df) == []

    def test_fires_one_finding_per_matching_row(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: x < 50, severity='warn',
            template='Sample {sample}: {value}',
        ))
        df = pd.DataFrame({
            'value': [10, 60, 30],
            'sample': ['a', 'b', 'c'],
        })
        out = evaluate_findings('m', df)
        assert len(out) == 2
        assert {f.sample for f in out} == {'a', 'c'}
        assert all(f.severity == 'warn' for f in out)

    def test_template_substitutes_row_columns(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: True, severity='warn',
            template='{sample} = {value} ({aa})',
        ))
        df = pd.DataFrame({'value': [42], 'sample': ['S1'], 'aa': ['Ala']})
        out = evaluate_findings('m', df)
        assert out[0].message == 'S1 = 42 (Ala)'

    def test_skips_nan_values(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: True, severity='warn',
            template='hit',
        ))
        df = pd.DataFrame({'value': [None, float('nan'), 1.0], 'sample': ['a', 'b', 'c']})
        out = evaluate_findings('m', df)
        # Only the third row triggers
        assert len(out) == 1
        assert out[0].sample == 'c'

    def test_highest_severity_wins_per_row(self, isolated_registry):
        # Two rules fire on the same row; keep the worst.
        register(FindingRule(
            anchor='m', predicate=lambda x: x < 50, severity='warn',
            template='warn',
        ))
        register(FindingRule(
            anchor='m', predicate=lambda x: x < 50, severity='fail',
            template='fail',
        ))
        df = pd.DataFrame({'value': [10], 'sample': ['a']})
        out = evaluate_findings('m', df)
        assert len(out) == 1
        assert out[0].severity == 'fail'
        assert out[0].message == 'fail'

    def test_template_failure_falls_back_to_raw(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: True, severity='warn',
            template='needs {missing_key}',
        ))
        df = pd.DataFrame({'value': [1], 'sample': ['a']})
        out = evaluate_findings('m', df)
        assert out[0].message == 'needs {missing_key}'  # raw template, no crash

    def test_handles_predicate_crash(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: 'not-a-bool' / 0, severity='warn',  # noqa: E501
            template='?',
        ))
        df = pd.DataFrame({'value': [1], 'sample': ['a']})
        assert evaluate_findings('m', df) == []

    def test_missing_value_column_returns_empty(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: True, severity='warn', template='t',
        ))
        df = pd.DataFrame({'other_col': [1]})
        assert evaluate_findings('m', df) == []

    def test_custom_value_column(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: x > 0, severity='warn',
            template='{sample}: {metric}',
        ))
        df = pd.DataFrame({'metric': [1, 0, 2], 'sample': ['a', 'b', 'c']})
        out = evaluate_findings('m', df, value_col='metric')
        assert len(out) == 2

    def test_type_error_when_not_dataframe(self, isolated_registry):
        with pytest.raises(TypeError):
            evaluate_findings('m', 'not a df')


# ---------------------------------------------------------------------------
# evaluate_value (single-value convenience)
# ---------------------------------------------------------------------------

class TestEvaluateValue:
    def test_single_value_fires(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: x < 50, severity='warn',
            template='Sample {sample}: {value:.1f}',
        ))
        f = evaluate_value('m', 30.0, sample='S1')
        assert f is not None
        assert f.severity == 'warn'
        assert 'S1' in f.message
        assert '30.0' in f.message

    def test_single_value_no_rule_returns_none(self, isolated_registry):
        assert evaluate_value('nonexistent', 10) is None

    def test_none_value_returns_none(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: True, severity='warn', template='t',
        ))
        assert evaluate_value('m', None) is None

    def test_picks_highest_severity(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: True, severity='warn', template='warn',
        ))
        register(FindingRule(
            anchor='m', predicate=lambda x: True, severity='fail', template='fail',
        ))
        f = evaluate_value('m', 1.0)
        assert f.severity == 'fail'

    def test_extra_kwargs_available_to_template(self, isolated_registry):
        register(FindingRule(
            anchor='m', predicate=lambda x: True, severity='warn',
            template='{value} on {chromosome}',
        ))
        f = evaluate_value('m', 5, extra={'chromosome': 'X'})
        assert 'X' in f.message


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

class TestAggregate:
    def test_groups_by_severity(self):
        fs = [
            _findings.Finding(anchor='a', severity='warn', message='m'),
            _findings.Finding(anchor='a', severity='warn', message='m'),
            _findings.Finding(anchor='a', severity='fail', message='m'),
            _findings.Finding(anchor='a', severity='pass', message='m'),
        ]
        d = aggregate_findings_by_severity(fs)
        assert len(d['warn']) == 2
        assert len(d['fail']) == 1
        assert len(d['pass']) == 1

    def test_empty_aggregate_has_all_keys(self):
        d = aggregate_findings_by_severity([])
        assert set(d.keys()) >= {'pass', 'warn', 'fail'}
        assert all(d[k] == [] for k in d)

    def test_summarize_findings_counts_and_worst(self):
        fs = [
            _findings.Finding(anchor='a', severity='warn', message='m'),
            _findings.Finding(anchor='a', severity='fail', message='m'),
        ]
        s = summarize_findings(fs)
        assert s == {'n_pass': 0, 'n_warn': 1, 'n_fail': 1, 'worst': 'fail'}

    def test_summarize_findings_clean_run(self):
        s = summarize_findings([])
        assert s == {'n_pass': 0, 'n_warn': 0, 'n_fail': 0, 'worst': 'pass'}
