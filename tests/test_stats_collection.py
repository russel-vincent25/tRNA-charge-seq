"""
Unit Tests for Chunked Stats Aggregation

Proves that the chunked read in :meth:`STATS_collection._collect_stats`
is output-identical to the previous whole-file ``pd.read_csv``. The
whole-file version materialised one row per aligned read (28 columns,
14 of them object-dtype strings), which OOM-killed workers at high
sequencing depth.

Run tests with:
    python -m pytest tests/test_stats_collection.py -v
"""

import bz2

import pandas as pd
import pytest

from src.stats_collection import STATS_collection


# ---------------------------------------------------------------------------
# Synthetic per-read stats CSV
# ---------------------------------------------------------------------------

# Low-cardinality values so that grouping keys genuinely collide both
# within and across chunk boundaries.
_ANNOTATIONS = ['Ecoli_tRNA-Ala-AGC-1-1', 'Ecoli_tRNA-Gly-GCC-1-1', 'Ecoli_tRNA-Lys-CTT-1-1']
_ANTICODONS = ['AGC', 'GCC', 'CTT']
_3P_NTS = ['CA', 'CC', 'GA']
# A mix of empty and non-empty so the row mask actually filters:
_3P_NON_TEMP = ['', '', '', 'GG', '']

_N_ROWS = 1500


def _synthetic_rows(n=_N_ROWS):
    """Yield ``n`` per-read stats records covering every 5p/3p_cover combination."""
    for i in range(n):
        ann_i = i % len(_ANNOTATIONS)
        yield {
            'readID': 'read_{}'.format(i),
            'common_seq': False,
            'sample_name_unique': 'PJ39-01-R1',
            'sample_name': 'PJ39-01',
            'replicate': 1,
            'barcode': 'BC01',
            'species': 'Ecoli',
            'tRNA_annotation': _ANNOTATIONS[ann_i],
            'align_score': 100 + (i % 7),
            'fmax_score': 0.95,
            'Ndeletions': i % 3,
            'Ninsertions': 0,
            'unique_annotation': True,
            'tRNA_annotation_len': 76,
            'align_5p_idx': 1,
            'align_3p_idx': 76,
            'align_5p_nt': 'G',
            'align_3p_nts': _3P_NTS[i % len(_3P_NTS)],
            'codon': 'GCT',
            'anticodon': _ANTICODONS[ann_i],
            'amino_acid': 'Ala',
            # Cycle all four 5p/3p_cover combinations so every fragment
            # class gets a non-zero count:
            '5p_cover': (i % 4) in (0, 1),
            '3p_cover': (i % 4) in (0, 2),
            '5p_non-temp': '',
            '3p_non-temp': _3P_NON_TEMP[i % len(_3P_NON_TEMP)],
            '5p_UMI': 'ACGTAC',
            '3p_BC': 'TTGG',
            'count': 1 + (i % 5),
        }


def _write_stats_csv(path, header, rows):
    """Write a bz2-compressed per-read stats CSV with the production header."""
    with bz2.open(path, 'wt') as fh:
        print(','.join(header), file=fh)
        for rec in rows:
            print(','.join(str(rec[col]) for col in header), file=fh)


# ---------------------------------------------------------------------------
# Reference implementation (the pre-chunking code path)
# ---------------------------------------------------------------------------

def _aggregate_whole_file(obj, stats_fnam):
    """
    The original whole-file implementation, verbatim, used as the oracle.

    Returns ``(fragment_counts, agg_df)``.
    """
    with bz2.open(stats_fnam, 'rt') as stats_fh:
        stat_df = pd.read_csv(stats_fh, keep_default_na=False, dtype=obj.stats_csv_header_td)

    ct = stat_df['count']
    _5p = stat_df['5p_cover'].astype(bool)
    _3p = stat_df['3p_cover'].astype(bool)
    fragment_counts = {
        'N_full_length': int(ct[_5p & _3p].sum()),
        'N_rt_dropoff': int(ct[~_5p & _3p].sum()),
        'N_5p_fragment': int(ct[_5p & ~_3p].sum()),
        'N_degraded': int(ct[~_5p & ~_3p].sum()),
        'N_total_aligned': int(ct.sum()),
    }

    row_mask = (stat_df['3p_cover']) & (stat_df['3p_non-temp'] == '')
    agg_df = stat_df[row_mask].groupby(obj.stats_agg_cols[:-1], as_index=False).agg({"count": "sum"})
    return fragment_counts, agg_df


def _sorted(df, obj):
    """Canonical ordering so frame comparison is order-insensitive."""
    return df.sort_values(obj.stats_agg_cols[:-1]).reset_index(drop=True)


def _read_agg_csv(path):
    """
    Read an aggregate CSV with no dtype coercion.

    Both sides of the comparison go through this so they are compared as
    the written artefact, which is what downstream stages actually consume.
    Note ``unique_annotation`` is typed ``str`` in ``stats_csv_header_td``
    but ``bool`` in ``stats_agg_cols_td``, so an in-memory oracle frame and
    a re-read output frame would otherwise differ on dtype alone.
    """
    return pd.read_csv(path, keep_default_na=False)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

class _StubbedSTATS(STATS_collection):
    """
    Drives the real ``_collect_stats`` against a canned per-read stats CSV.

    ``_read_non_common`` normally streams the SWalign JSON and the trimmed
    FASTQ; here it just emits the synthetic records, so the test exercises
    the production aggregation path without needing pipeline inputs.
    """

    rows_to_emit = []

    def _read_non_common(self, row, stats_fh):
        for rec in self.rows_to_emit:
            print(','.join(str(rec[col]) for col in self.stats_csv_header), file=stats_fh)


@pytest.fixture
def stats_obj(tmp_path):
    """A STATS_collection wired to a temporary project directory."""
    data_dir = tmp_path / 'data'
    (data_dir / 'align').mkdir(parents=True)
    dir_dict = {
        'NBdir': str(tmp_path),
        'data_dir': 'data',
        'align_dir': 'align',
        'UMI_dir': 'UMI',
        'stats_dir': 'stats',
    }
    sample_df = pd.DataFrame([{
        'sample_name_unique': 'PJ39-01-R1',
        'sample_name': 'PJ39-01',
        'replicate': 1,
        'barcode': 'BC01',
        'species': 'Ecoli',
    }])
    obj = _StubbedSTATS(dir_dict, tRNA_data=dict(), sample_df=sample_df,
                        check_exists=False, overwrite_dir=True)
    obj.verbose = False  # normally set by run_parallel()
    return obj


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestChunkedAggregation:
    """The chunked read must be byte-for-byte equivalent to the whole-file read."""

    # Deliberately includes sizes that do NOT divide 1500 evenly, one larger
    # than the row count (single chunk), and 1 (every row its own chunk).
    @pytest.mark.parametrize('chunksize', [1, 7, 100, 333, 749, 1499, 1500, 5000])
    def test_matches_whole_file_read(self, stats_obj, chunksize):
        stats_obj.rows_to_emit = list(_synthetic_rows())
        stats_obj.stats_chunksize = chunksize

        row = stats_obj.sample_df.iloc[0]
        result = stats_obj._collect_stats(0, row)

        stats_fnam = '{}/{}_stats.csv.bz2'.format(stats_obj.stats_dir_abs, row['sample_name_unique'])
        expected_counts, expected_agg = _aggregate_whole_file(stats_obj, stats_fnam)

        assert result['fragment_counts'] == expected_counts

        expected_path = '{}/expected_aggregate.csv'.format(stats_obj.stats_dir_abs)
        expected_agg.to_csv(expected_path, header=True, index=False)

        pd.testing.assert_frame_equal(_sorted(_read_agg_csv(result['stats_agg_path']), stats_obj),
                                      _sorted(_read_agg_csv(expected_path), stats_obj))

    def test_fragment_counts_are_non_trivial(self, stats_obj):
        """Guard the oracle: every fragment class must be exercised."""
        stats_obj.rows_to_emit = list(_synthetic_rows())
        stats_obj.stats_chunksize = 250

        row = stats_obj.sample_df.iloc[0]
        counts = stats_obj._collect_stats(0, row)['fragment_counts']

        for key in ('N_full_length', 'N_rt_dropoff', 'N_5p_fragment', 'N_degraded'):
            assert counts[key] > 0, key
        assert counts['N_total_aligned'] == (counts['N_full_length'] + counts['N_rt_dropoff']
                                             + counts['N_5p_fragment'] + counts['N_degraded'])

    def test_row_mask_actually_filters(self, stats_obj):
        """The aggregate must drop 3p_non-temp rows, not just pass everything through."""
        stats_obj.rows_to_emit = list(_synthetic_rows())
        stats_obj.stats_chunksize = 250

        row = stats_obj.sample_df.iloc[0]
        result = stats_obj._collect_stats(0, row)
        agg = _read_agg_csv(result['stats_agg_path'])

        assert agg['count'].sum() < result['fragment_counts']['N_total_aligned']

    @pytest.mark.parametrize('chunksize', [1, 500000])
    def test_header_only_input(self, stats_obj, chunksize):
        """A stats CSV with no data rows still yields a correctly headed empty aggregate."""
        stats_obj.rows_to_emit = []
        stats_obj.stats_chunksize = chunksize

        row = stats_obj.sample_df.iloc[0]
        result = stats_obj._collect_stats(0, row)

        assert result['fragment_counts'] == {
            'N_full_length': 0,
            'N_rt_dropoff': 0,
            'N_5p_fragment': 0,
            'N_degraded': 0,
            'N_total_aligned': 0,
        }

        agg = _read_agg_csv(result['stats_agg_path'])
        assert len(agg) == 0
        assert list(agg.columns) == stats_obj.stats_agg_cols

    def test_default_chunksize(self, stats_obj):
        assert stats_obj.stats_chunksize == 500000
