"""QC reporting module for tRNA-charge-seq pipeline."""

from trnaseq.qc.report import QCReportGenerator

try:
    from trnaseq.qc.modification_report import ModificationReportGenerator
except ImportError:
    ModificationReportGenerator = None

__all__ = ['QCReportGenerator', 'ModificationReportGenerator']
