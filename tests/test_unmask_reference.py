"""Tests for restoring N-masked reference bases before rate calculation."""

import pandas as pd
import pytest

from trnaseq.modifications.positional import unmask_reference
from trnaseq.modifications.rt_signatures import RTSignatureAnalyzer


def _fasta(path, records):
    path.write_text(''.join(f'>{n}\n{s}\n' for n, s in records.items()))
    return path


def test_restores_masked_bases(tmp_path):
    ref = {'t1': {'seq': 'ACNGTN', 'seq_len': 6}, 't2': {'seq': 'ACGT', 'seq_len': 4}}
    fa = _fasta(tmp_path / 'u.fa', {'t1': 'ACAGTC', 't2': 'ACGT'})
    assert unmask_reference(ref, fa) == 2
    assert ref['t1']['seq'] == 'ACAGTC' and ref['t2']['seq'] == 'ACGT'


def test_rejects_a_fasta_that_is_not_the_counterpart(tmp_path):
    ref = {'t1': {'seq': 'ACNGTN', 'seq_len': 6}}
    with pytest.raises(ValueError):
        unmask_reference(ref, _fasta(tmp_path / 'u.fa', {'t1': 'TCAGTC'}))
    with pytest.raises(ValueError):
        unmask_reference(ref, _fasta(tmp_path / 'v.fa', {'t1': 'ACAGT'}))


def test_missing_reference_warns_and_stays_masked(tmp_path):
    ref = {'t1': {'seq': 'ACNG', 'seq_len': 4}}
    with pytest.warns(UserWarning):
        assert unmask_reference(ref, _fasta(tmp_path / 'u.fa', {'other': 'ACAG'})) == 0
    assert ref['t1']['seq'] == 'ACNG'


def test_substitutions_not_measured_against_n():
    analyzer = RTSignatureAnalyzer(min_coverage=10, verbose=False)
    pscm = pd.DataFrame({'A': [100, 0], 'C': [0, 0], 'G': [0, 0], 'T': [0, 100],
                         'U': [0, 0], 'N': [0, 0], '-': [0, 0]})
    masked = analyzer.calculate_mismatch_rates(pscm, 'AN')
    restored = analyzer.calculate_mismatch_rates(pscm, 'AA')
    assert masked['mismatch_rate'].tolist() == [0.0, 0.0]   # N: not measured
    assert restored['mismatch_rate'].tolist() == [0.0, 1.0]  # real base: A->T measured
