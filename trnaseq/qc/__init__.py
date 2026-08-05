"""QC reporting module for tRNA-charge-seq pipeline."""

from trnaseq.qc._common import (
    ReportContext,
    render_html_shell,
    render_panel,
    render_skipped_panel,
    fig_to_div,
    export_panel_data,
    resolve_trna_col,
    host_filter,
    categorical_palette,
    OKABE_ITO,
    GLASBEY_HEX,
    STYLE_CSS,
)
from trnaseq.qc._findings import (
    FindingRule,
    Finding,
    FINDINGS,
    register as register_finding,
    rules_for as findings_for,
    evaluate_findings,
    evaluate_value,
    aggregate_findings_by_severity,
    worst_severity,
    summarize_findings,
)
from trnaseq.qc.report import QCReportGenerator

try:
    from trnaseq.qc.modification_report import ModificationReportGenerator
except ImportError:
    ModificationReportGenerator = None

try:
    from trnaseq.qc.charge_report import ChargeReportGenerator
except ImportError:
    ChargeReportGenerator = None

try:
    from trnaseq.qc.fragment_report import FragmentReportGenerator
except ImportError:
    FragmentReportGenerator = None

try:
    from trnaseq.qc.abundance_report import AbundanceReportGenerator
except ImportError:
    AbundanceReportGenerator = None

__all__ = [
    # Generators
    'QCReportGenerator',
    'ModificationReportGenerator',
    'ChargeReportGenerator',
    'FragmentReportGenerator',
    'AbundanceReportGenerator',
    # Shared infrastructure (T0)
    'ReportContext',
    'render_html_shell',
    'render_panel',
    'render_skipped_panel',
    'fig_to_div',
    'export_panel_data',
    'resolve_trna_col',
    'host_filter',
    'categorical_palette',
    'OKABE_ITO',
    'GLASBEY_HEX',
    'STYLE_CSS',
    # Findings registry (T0)
    'FindingRule',
    'Finding',
    'FINDINGS',
    'register_finding',
    'findings_for',
    'evaluate_findings',
    'evaluate_value',
    'aggregate_findings_by_severity',
    'worst_severity',
    'summarize_findings',
]
