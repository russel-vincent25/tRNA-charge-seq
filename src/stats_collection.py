import os, shutil, bz2, warnings, json, gc
import json_stream
from Bio import SeqIO
import pandas as pd
import numpy as np
from mpire import WorkerPool



class STATS_collection:
    '''
    Class to collect statistics from the
    alignment results.
    Keyword arguments:
    common_seqs -- bzip2 compressed fasta file of commonly observed sequences to avoid duplicated alignments (default None)
    ignore_common_count -- Ignore common count even if X_common-seq-obs.json filename exists (default False)
    check_exists -- Check if required files exist before starting. If set to False, this will ignore checking for common count (default True)
    overwrite_dir -- Overwrite old stats folder if any exists (default False)
    check_exists -- Check if input files exist (default True)
    reads_SW_sorted -- Assume reads and SW results are sorted in the same order. This gives a massive memory saving. (default True)
    from_UMIdir -- Is the input data from a folder made with the UMI_trim class? (default True)
    '''
    def __init__(self, dir_dict, tRNA_data, sample_df, common_seqs=None, \
                 ignore_common_count=False, check_exists=True, overwrite_dir=False, \
                 reads_SW_sorted=True, from_UMIdir=True):
        self.stats_csv_header = ['readID', 'common_seq', 'sample_name_unique', \
                                 'sample_name', 'replicate', 'barcode', 'species', 'tRNA_annotation', \
                                 'align_score', 'fmax_score', 'Ndeletions', 'Ninsertions', \
                                 'unique_annotation', 'tRNA_annotation_len', \
                                 'align_5p_idx', 'align_3p_idx', 'align_5p_nt', 'align_3p_nts', \
                                 'codon', 'anticodon', 'amino_acid', '5p_cover', '3p_cover', \
                                 '5p_non-temp', '3p_non-temp', '5p_UMI', '3p_BC', 'count']
        self.stats_csv_header_type = [str, bool, str, str, int, str, str, str, int, float, \
                                      int, int, str, int, int, int, str, str, str, str, str, \
                                      bool, bool, str, str, str, str, int]
        self.stats_csv_header_td = {nam:tp for nam, tp in zip(self.stats_csv_header, self.stats_csv_header_type)}
        # Here: could add number of gaps, or maybe a boolean, indicating if the faction or align score to max is above a threshold
        self.stats_agg_cols = ['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'species', \
                               'tRNA_annotation', 'tRNA_annotation_len', 'unique_annotation', \
                               '5p_cover', 'align_3p_nts', 'codon', 'anticodon', 'amino_acid', \
                               'count']
        self.stats_agg_cols_type = [str, str, int, str, str, str, int, bool, \
                                    bool, str, str, str, str, int]
        self.stats_agg_cols_td = {nam:tp for nam, tp in zip(self.stats_agg_cols, self.stats_agg_cols_type)}

        # Input:
        self.tRNA_data, self.sample_df = tRNA_data, sample_df
        self.dir_dict = dir_dict
        self.common_seqs_fnam = common_seqs
        self.common_seqs_info = dict() # Name to sequence for common sequence
        self.reads_SW_sorted = reads_SW_sorted
        self.from_UMIdir = from_UMIdir

        # Check files exists before starting:
        self.align_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['align_dir'])
        if self.from_UMIdir:
            self.UMI_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['UMI_dir'])
        if check_exists:
            for _, row in self.sample_df.iterrows():
                SWres_fnam = '{}/{}_SWalign.json.bz2'.format(self.align_dir_abs, row['sample_name_unique'])
                assert(os.path.exists(SWres_fnam))
                if self.from_UMIdir:
                    trimmed_fn = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
                    assert(os.path.exists(trimmed_fn))
                else:
                    # File paths are specified in sample_df
                    for _, row in self.sample_df.iterrows():
                        # Absolute path:
                        if row['path'][0] == '/':
                            fpath = row['path']
                        # Relative path:
                        else:
                            fpath = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], row['path'])
                    assert(os.path.exists(fpath))
                common_obs_fn = '{}/{}_common-seq-obs.json'.format(self.align_dir_abs, row['sample_name_unique'])
                if not self.common_seqs_fnam is None:
                    assert(os.path.exists(common_obs_fn))
                elif self.common_seqs_fnam is None and ignore_common_count is False and os.path.exists(common_obs_fn):
                    raise Exception('Found common sequence counts for {}\n'
                                    'See: {}\n'
                                    'But common sequences are not specified. Either specify common sequences '
                                    'or explicitly set ignore_common_count=True'.format(row['sample_name_unique'], common_obs_fn))

        # Make output folder:
        self._make_dir(overwrite=overwrite_dir)

    def _make_dir(self, overwrite):
        # Create folder for files:
        self.stats_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['stats_dir'])
        try:
            os.mkdir(self.stats_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.stats_dir_abs)
                os.mkdir(self.stats_dir_abs)
            else:
                print('Using existing folder because overwrite set to false: {}'.format(self.stats_dir_abs))

    def run_parallel(self, n_jobs=4, verbose=True, \
                     load_previous=False):
        '''
        Submit the input files for stats collection.
        Keyword arguments:
        n_jobs -- Number of subprocesses started in parallel (default 4)
        load_previous -- Attempt to load results from a previous stats collection by looking up ALL_stats_aggregate.csv (default False)
        verbose -- Verbose printing collection progress (default True)
        '''
        if load_previous:
            stats_agg_fnam = '{}/ALL_stats_aggregate.csv'.format(self.stats_dir_abs)
            self.concat_df = pd.read_csv(stats_agg_fnam, keep_default_na=False, dtype=self.stats_agg_cols_td)
            print('Loaded results from previous run... Not running stats collection.')
            return(self.concat_df)
        elif not self.common_seqs_fnam is None:
            self._load_common_seqs(verbose)

        self.verbose = verbose
        if self.verbose:
            print('Collecting stats from:', end='')
        # Run parallel:
        data = list(self.sample_df.iterrows())
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self._collect_stats, data)
        self._concat_stats(results)
        return(self.concat_df)

    '''
    def run_serial(self, verbose=True):
        self.verbose = verbose
        if self.verbose:
            print('Collecting stats from:', end='')
        results = [self._collect_stats(index, row) for index, row in self.sample_df.iterrows()]
        self._concat_stats(results)
        return(self.concat_df)
    '''

    def _load_common_seqs(self, verbose):
        # Make name to sequence dictionary for common sequences.
        # We can only allow one species if using common sequences.
        # Multiple species would require running the alignment on common sequences
        # several times, defeating the purpose, but also making the code much
        # more complicated.
        sp_set = set(self.sample_df['species'].values)
        if len(sp_set) > 1:
            raise Exception('Only one species allowed in sample sheet when using common sequences.')
        self.common_seqs_sp = list(sp_set)[0]

        if verbose:
            print('Using common sequences...')
        assert(self.common_seqs_fnam[-4:] == '.bz2')
        with bz2.open(self.common_seqs_fnam, "rt") as input_fh:
            for ridx, record in enumerate(SeqIO.parse(input_fh, "fasta")):
                assert(ridx == int(record.id))
                self.common_seqs_info[record.id] = str(record.seq)

    def _collect_stats(self, index, row):
        if self.verbose:
            print('  {}'.format(row['sample_name_unique']), end='')

        # Write stats to file:
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        stats_agg_fnam = '{}/{}_stats_aggregate.csv'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'wt') as stats_fh:
            # Print header to stats CSV file:
            print(','.join(self.stats_csv_header), file=stats_fh)
            self._read_non_common(row, stats_fh)
            if not self.common_seqs_fnam is None:
                self._read_common(row, stats_fh)

        # Filter data and aggregate to count charged/uncharged tRNAs
        # Read stats from stats CSV file:
        with bz2.open(stats_fnam, 'rt') as stats_fh:
            # Use "keep_default_na=False" to read an empty string
            # as an empty string and not as NaN.
            stat_df = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)

        # Compute fragment counts BEFORE the 3p_cover filter
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

        # Aggregate dataframe and write as CSV file:
        # Here: also filter sequences with long 5p_non-temp sequences (these are likely template switch products)
        row_mask = (stat_df['3p_cover']) & (stat_df['3p_non-temp'] == '')
        agg_df = stat_df[row_mask].groupby(self.stats_agg_cols[:-1], as_index=False).agg({"count": "sum"})
        agg_df.to_csv(stats_agg_fnam, header=True, index=False)

        return({
            'stats_agg_path': stats_agg_fnam,
            'fragment_counts': fragment_counts,
            'sample_name_unique': row['sample_name_unique'],
        })

    def _read_non_common(self, row, stats_fh):
        # Read fastq files must be read to
        # extract UMI and 5/3p non-template bases:
        # File from UMI dir or path specified in sample_df:
        if self.from_UMIdir:
            trimmed_fn = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
        else:
            # Absolute path:
            if row['path'][0] == '/':
                trimmed_fn = row['path']
            # Relative path:
            else:
                trimmed_fn = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], row['path'])
        if self.reads_SW_sorted:
            reads_fh = bz2.open(trimmed_fn, 'rt')
            read_info_iter = SeqIO.parse(reads_fh, "fastq")
            readID = ''
        # If reads and SW results are not sorted the same way,
        # the read information must be read into memory:
        else:
            read_info = dict()
            with bz2.open(trimmed_fn, 'rt') as fh_bz:
                for read in SeqIO.parse(fh_bz, "fastq"):
                    # The last two strings are the adapter sequence and the UMI:
                    try:
                        _3p_bc, _5p_umi = read.description.split()[-1].split(':')[-2:]
                    except:
                        _3p_bc, _5p_umi = '', ''
                    seq = str(read.seq)
                    readID = str(read.id)
                    read_info[readID] = (_5p_umi, seq)

        # Open the alignment results:
        SWres_fnam = '{}/{}_SWalign.json.bz2'.format(self.align_dir_abs, row['sample_name_unique'])
        with bz2.open(SWres_fnam, 'rt', encoding="utf-8") as SWres_fh:
            # Parse JSON data as a stream (saves memory),
            # i.e. as a transient dict-like object
            SWres = json_stream.load(SWres_fh)
            # Loop through each read in the alignment results:
            for SWreadID, align_dict in SWres.persistent().items():
                common_seq = False
                # Skip reads that were not aligned:
                if not align_dict['aligned']:
                    continue

                # Extract read info:
                try:
                    if self.reads_SW_sorted:
                        while readID != SWreadID:
                            read = next(read_info_iter)
                            try:
                                _3p_bc, _5p_umi = read.description.split()[-1].split(':')[-2:]
                            except:
                                _3p_bc, _5p_umi = '', ''
                            read_seq = str(read.seq)
                            readID = str(read.id)
                    else:
                        _5p_umi, read_seq = read_info.pop(SWreadID)
                except:
                    raise Exception('Read ID ({}) not found among read sequences. Did any of the fastq headers change such that there is a mismatch between headers in the alignment json and those in the reads?'.format(SWreadID))

                # Collect all the information:
                sample_name_unique = row['sample_name_unique']
                species = row['species']
                sample_name = row_exist_or_none(row, 'sample_name')
                replicate = row_exist_or_none(row, 'replicate')
                barcode = row_exist_or_none(row, 'barcode')
                tRNA_annotation = align_dict['name']
                tRNA_annotation_first = tRNA_annotation.split('@')[0]
                align_score = align_dict['score']
                fmax_score = align_dict['Fmax_score']
                Ndel = align_dict['Ndel']
                Nins = align_dict['Nins']
                unique_annotation = '@' not in tRNA_annotation
                tRNA_annotation_len = self.tRNA_data[tRNA_annotation_first]['len']
                align_5p_idx, align_3p_idx = align_dict['dpos']
                align_5p_nt = align_dict['qseq'][0]
                qseq = align_dict['qseq']
                align_3p_nts = qseq[-2:] if len(qseq) >= 2 else qseq

                # Move index for reads with 3' cleaved A:
                if align_3p_idx == (tRNA_annotation_len - 1) and align_3p_nts[-1] == 'C':
                    align_3p_idx += 1
                codon = self.tRNA_data[tRNA_annotation_first]['codon']
                anticodon = self.tRNA_data[tRNA_annotation_first]['anticodon']
                amino_acid = self.tRNA_data[tRNA_annotation_first]['amino_acid']
                _5p_cover = align_5p_idx == 1
                _3p_cover = align_3p_idx == tRNA_annotation_len
                qpos = align_dict['qpos']
                _5p_non_temp = read_seq[0:(qpos[0]-1)]
                _3p_non_temp = read_seq[qpos[1]:]
                _3p_bc = row_exist_or_none(row, 'barcode_seq')
                # For "non-common" sequences multiple reads
                # have not been collapsed:
                count = 1

                # Print line to output csv file:
                line_lst = [readID, common_seq, sample_name_unique, sample_name, replicate, \
                            barcode, species, tRNA_annotation, align_score, fmax_score, \
                            Ndel, Nins, unique_annotation, \
                            tRNA_annotation_len, align_5p_idx, align_3p_idx, align_5p_nt, \
                            align_3p_nts, codon, anticodon, amino_acid, _5p_cover, _3p_cover, \
                            _5p_non_temp, _3p_non_temp, _5p_umi, _3p_bc, count]
                csv_line = ','.join(map(str, line_lst))
                print(csv_line, file=stats_fh)

        reads_fh.close()

    def _read_common(self, row, stats_fh):
        # Read common sequences observations for this sample:
        common_obs_fn = '{}/{}_common-seq-obs.json'.format(self.align_dir_abs, row['sample_name_unique'])
        with open(common_obs_fn, 'r') as fh_in:
            common_obs_raw = json.load(fh_in)
        # Handle both formats: raw list or dict with 'common_obs' key
        if isinstance(common_obs_raw, dict) and 'common_obs' in common_obs_raw:
            common_obs = common_obs_raw['common_obs']
        else:
            common_obs = common_obs_raw

        # Open the alignment results:
        SWres_fnam = '{}/{}_SWalign.json.bz2'.format(self.align_dir_abs, 'common-seqs')
        with bz2.open(SWres_fnam, 'rt', encoding="utf-8") as SWres_fh:
            # Parse JSON data as a stream (saves memory),
            # i.e. as a transient dict-like object
            SWres = json_stream.load(SWres_fh)
            # Loop through each read in the alignment results:
            for readID, align_dict in SWres.persistent().items():
                common_seq = True
                readID_int = int(readID)
                # Skip reads that were not aligned:
                if not align_dict['aligned']:
                    continue
                # And skip if this sample did not have this common sequence:
                elif common_obs[readID_int] == 0:
                    continue

                # Collect all the information:
                sample_name_unique = row['sample_name_unique']
                species = row['species']
                sample_name = row_exist_or_none(row, 'sample_name')
                replicate = row_exist_or_none(row, 'replicate')
                barcode = row_exist_or_none(row, 'barcode')
                tRNA_annotation = align_dict['name']
                tRNA_annotation_first = tRNA_annotation.split('@')[0]
                align_score = align_dict['score']
                fmax_score = align_dict['Fmax_score']
                Ndel = align_dict['Ndel']
                Nins = align_dict['Nins']
                unique_annotation = '@' not in tRNA_annotation
                tRNA_annotation_len = self.tRNA_data[tRNA_annotation_first]['len']
                align_5p_idx, align_3p_idx = align_dict['dpos']
                align_5p_nt = align_dict['qseq'][0]
                qseq = align_dict['qseq']
                align_3p_nts = qseq[-2:] if len(qseq) >= 2 else qseq

                # Move index for reads with 3' cleaved A:
                if align_3p_idx == (tRNA_annotation_len - 1) and align_3p_nts[-1] == 'C':
                    align_3p_idx += 1
                codon = self.tRNA_data[tRNA_annotation_first]['codon']
                anticodon = self.tRNA_data[tRNA_annotation_first]['anticodon']
                amino_acid = self.tRNA_data[tRNA_annotation_first]['amino_acid']
                _5p_cover = align_5p_idx == 1
                _3p_cover = align_3p_idx == tRNA_annotation_len

                # Extract non-template bases from common reads:
                try:
                    seq = self.common_seqs_info[readID]
                except KeyError:
                    raise Exception('Read ID ({}) not found among sequences. '
                                    'Did any of the fastq headers change such that there is '
                                    'a mismatch between headers in the alignment json and those '
                                    'in the reads?'.format(readID))
                qpos = align_dict['qpos']
                _5p_non_temp = seq[0:(qpos[0]-1)]
                _3p_non_temp = seq[qpos[1]:]
                _5p_umi = ''  # UMI information is lost when using common sequences
                _3p_bc = row_exist_or_none(row, 'barcode_seq')
                # For common sequences the add the read count:
                count = int(common_obs[readID_int])

                # Print line to output csv file:
                line_lst = [readID, common_seq, sample_name_unique, sample_name, replicate, \
                            barcode, species, tRNA_annotation, align_score, fmax_score, \
                            Ndel, Nins, unique_annotation, \
                            tRNA_annotation_len, align_5p_idx, align_3p_idx, align_5p_nt, \
                            align_3p_nts, codon, anticodon, amino_acid, _5p_cover, _3p_cover, \
                            _5p_non_temp, _3p_non_temp, _5p_umi, _3p_bc, count]
                csv_line = ','.join(map(str, line_lst))
                print(csv_line, file=stats_fh)

    def _concat_stats(self, results):
        # Backward compat: results may be list of strings (legacy) or list of dicts (new)
        if results and isinstance(results[0], str):
            csv_paths = results
            self.fragment_counts = {}
        else:
            csv_paths = [r['stats_agg_path'] for r in results]
            self.fragment_counts = {r['sample_name_unique']: r['fragment_counts'] for r in results}

        # Concatenate all the aggregated stats csv files:
        stats_agg_fnam = '{}/ALL_stats_aggregate.csv'.format(self.stats_dir_abs)
        with open(stats_agg_fnam, 'w') as fh_out:
            # Grap header:
            with open(csv_paths[0], 'r') as fh_in:
                print(fh_in.readline(), file=fh_out, end='')
            for path in csv_paths:
                with open(path, 'r') as fh_in:
                    next(fh_in) # burn the header
                    for line in fh_in:
                        print(line, file=fh_out, end='')

        # Read and store as dataframe:
        self.concat_df = pd.read_csv(stats_agg_fnam, keep_default_na=False, dtype=self.stats_agg_cols_td)


def row_exist_or_none(row, col):
    try:
        return(row[col])
    except:
        return(None)
