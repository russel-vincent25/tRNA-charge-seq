"""QC reporting module for tRNA-charge-seq pipeline."""

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
    'QCReportGenerator',
    'ModificationReportGenerator',
    'ChargeReportGenerator',
    'FragmentReportGenerator',
    'AbundanceReportGenerator',
]
