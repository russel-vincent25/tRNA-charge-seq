import os, shutil, bz2, resource, re
from collections import Counter
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import numpy as np
from mpire import WorkerPool
import jellyfish



class AR_merge:
    '''
    Class to merge the paired end reads using AdapterRemoval.

    Keyword arguments:
    AR_threads -- Threads specified to AdapterRemoval (default 4)
    overwrite_dir -- Overwrite old merge folder if any exists (default False)
    check_input -- Check if input files exist (default True)
    '''
    def __init__(self, dir_dict, inp_file_df, MIN_READ_LEN, \
                 AR_threads=4, overwrite_dir=False, check_input=True):
        # Input:
        self.inp_file_df, self.MIN_READ_LEN = inp_file_df, MIN_READ_LEN
        self.dir_dict = dir_dict
        # AdapterRomoval input:
        self.adapter1_tmp = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC<P7_index>ATCTCGTATGCCGTCTTCTGCTTG'
        self.adapter2_tmp = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT<P5_index>GTGTAGATCTCGGTGGTCGCCGTATCATT'
        self.AR_cmd_tmp = ["AdapterRemoval", "--bzip2", "--preserve5p", "--collapse", "--minalignmentlength", "10", "--threads", str(AR_threads)]
        self.AR_overwrite = True

        self.seq_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['seq_dir'])
        # Check files exists before starting:
        if check_input:
            for _, row in self.inp_file_df.iterrows():
                fnam_mate1 = '{}/{}'.format(self.seq_dir_abs, row['fastq_mate1_filename'])
                fnam_mate2 = '{}/{}'.format(self.seq_dir_abs, row['fastq_mate2_filename'])
                assert(os.path.exists(fnam_mate1))
                assert(os.path.exists(fnam_mate2))

        # Make output folder:
        self._make_dir(overwrite_dir)

    def _make_dir(self, overwrite=True):
        # Create folder for files:
        self.AdapterRemoval_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['AdapterRemoval_dir'])
        try:
            os.mkdir(self.AdapterRemoval_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.AdapterRemoval_dir_abs)
                os.mkdir(self.AdapterRemoval_dir_abs)
            else:
                print('Using existing folder because overwrite set to false: {}'.format(self.AdapterRemoval_dir_abs))

    def run_parallel(self, n_jobs=4, overwrite=True):
        '''
        Submit the mate pair files for merging.

        Keyword arguments:
        n_jobs -- Number of subprocesses of AdapterRemoval started in parallel (default 4)
        overwrite -- Overwrite files of previous run. If false, skipping mate pair files with merged files existing (default True)
        '''
        self.AR_overwrite = overwrite
        os.chdir(self.AdapterRemoval_dir_abs)
        try:
            data = list(self.inp_file_df.iterrows())
            with WorkerPool(n_jobs=n_jobs) as pool:
                results = pool.map(self._start_AR, data)
            self._collect_stats()
            os.chdir(self.dir_dict['NBdir'])
            return(self.inp_file_df)
        except Exception as err:
            os.chdir(self.dir_dict['NBdir'])
            raise err

    def _start_AR(self, index, row):
        AR_cmd = self.AR_cmd_tmp.copy()
        basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
        # Check if output exists, and skip if not overwrite:
        merged_fastq_fn = '{}.collapsed.bz2'.format(basename)
        log_fn = '{}_logfile.txt'.format(basename)
        settings_fn = '{}.settings'.format(basename)
        if not self.AR_overwrite and os.path.isfile(merged_fastq_fn) and os.path.isfile(log_fn) and os.path.isfile(settings_fn):
            return(1)
        adapter1 = self.adapter1_tmp.replace('<P7_index>', row['P7_index_seq'])
        adapter2 = self.adapter2_tmp.replace('<P5_index>', row['P5_index_seq'])
        fnam_mate1 = '{}/{}'.format(self.seq_dir_abs, row['fastq_mate1_filename'])
        fnam_mate2 = '{}/{}'.format(self.seq_dir_abs, row['fastq_mate2_filename'])
        AR_cmd.extend(['--minlength', str(self.MIN_READ_LEN)])
        AR_cmd.extend(['--adapter1', adapter1])
        AR_cmd.extend(['--adapter2', adapter2])
        AR_cmd.extend(['--basename', basename])
        AR_cmd.extend(['--file1', '{}'.format(fnam_mate1)])
        AR_cmd.extend(['--file2', '{}'.format(fnam_mate2)])

        # Run AdapterRemoval, collect log in "file":
        with Popen(AR_cmd, stdout=PIPE, stderr=STDOUT) as p, open(log_fn, 'a') as file:
            file.write('Starting subprocess with command:')
            file.write(str(AR_cmd))
            file.write('\n')
            for line in p.stdout:
                file.write(line.decode('utf-8'))
            file.write('\n****** DONE ******\n\n\n')
        return(1)

    def _collect_stats(self):
        N_pairs = list()
        N_merged = list()
        for _, row in self.inp_file_df.iterrows():
            basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
            with open('{}.settings'.format(basename), 'r') as fh:
                for line in fh:
                    if 'Total number of read pairs:' in line:
                        N_pairs.append(int(line.split(':')[1][1:]))
                    if 'Number of full-length collapsed pairs:' in line:
                        N_merged.append(int(line.split(':')[1][1:]))
        # Write stats:
        self.inp_file_df['N_pairs'] = N_pairs
        self.inp_file_df['N_merged'] = N_merged
        self.inp_file_df['percent_successfully_merged'] = self.inp_file_df['N_merged'].values / self.inp_file_df['N_pairs'].values *100
        self.inp_file_df.to_excel('merge_stats.xlsx')

    def qc_probe(self, barcodes, umi_mode='anchored', anchor_seq=None,
                 max_stagger=3, anchor_max_dist=1, umi_len=10,
                 umi_end={'T', 'C'}, n_reads=10000):
        """
        Lightweight QC probe run on merged reads to check UMI and barcode
        detectability *before* committing to full BC_split / UMI_trim stages.

        Parameters:
            barcodes: list of barcode sequences to look for at 3' end
            umi_mode: 'anchored' or 'pyrimidine'
            anchor_seq: constant anchor sequence (required when umi_mode='anchored')
            max_stagger: max random nt before anchor (anchored mode)
            anchor_max_dist: hamming tolerance for anchor match
            umi_len: UMI length (pyrimidine mode)
            umi_end: set of valid last-nt for UMI (pyrimidine mode)
            n_reads: number of reads to subsample per file pair

        Returns:
            dict mapping basename -> {n_probed, pct_umi, pct_bc, pct_both,
                                      stagger_counts (anchored mode only)}
        """
        results = {}
        anchor_len = len(anchor_seq) if anchor_seq else 0

        for _, row in self.inp_file_df.iterrows():
            basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
            merged_fn = '{}/{}.collapsed.bz2'.format(self.AdapterRemoval_dir_abs, basename)
            if not os.path.exists(merged_fn):
                continue

            n_probed = 0
            n_umi = 0
            n_bc = 0
            n_both = 0
            stagger_counts = Counter()

            with bz2.open(merged_fn, "rt") as fh:
                for title, seq, qual in FastqGeneralIterator(fh):
                    if n_probed >= n_reads:
                        break
                    n_probed += 1

                    # --- 5' UMI check ---
                    has_umi = False
                    if umi_mode == 'anchored' and anchor_seq:
                        for stagger in range(max_stagger + 1):
                            if stagger + anchor_len > len(seq):
                                break
                            candidate = seq[stagger:stagger + anchor_len]
                            if jellyfish.hamming_distance(candidate, anchor_seq) <= anchor_max_dist:
                                has_umi = True
                                stagger_counts[stagger] += 1
                                break
                    elif umi_mode == 'pyrimidine':
                        if len(seq) >= umi_len and seq[umi_len - 1] in umi_end:
                            has_umi = True

                    # --- 3' barcode check ---
                    has_bc = False
                    for bc in barcodes:
                        bc_len = len(bc)
                        if bc_len > len(seq):
                            continue
                        tail = seq[-bc_len:]
                        if 'N' in bc:
                            mismatches = sum(a != b for a, b in zip(tail, bc) if b != 'N')
                            if mismatches <= 1:
                                has_bc = True
                                break
                        else:
                            if jellyfish.hamming_distance(tail, bc) <= 1:
                                has_bc = True
                                break

                    if has_umi:
                        n_umi += 1
                    if has_bc:
                        n_bc += 1
                    if has_umi and has_bc:
                        n_both += 1

            result = {
                'n_probed': n_probed,
                'pct_umi': round(n_umi / n_probed * 100, 1) if n_probed else 0,
                'pct_bc': round(n_bc / n_probed * 100, 1) if n_probed else 0,
                'pct_both': round(n_both / n_probed * 100, 1) if n_probed else 0,
            }
            if umi_mode == 'anchored':
                result['stagger_counts'] = dict(stagger_counts)
            results[basename] = result

        return results



class BC_split:
    '''
    Class to split the merged fastq files into several files based on barcodes.

    Keyword arguments:
    max_dist -- Maximum hamming distance between a read and a barcode to assign barcode identify to the read (default 1)
    overwrite_dir -- Overwrite old barcode split folder if any exists (default False)
    '''
    def __init__(self, dir_dict, sample_df, inp_file_df, overwrite_dir=False, \
                 max_dist=1):
        # Input:
        self.sample_df, self.inp_file_df = sample_df, inp_file_df
        self.dir_dict = dir_dict
        self.max_dist = max_dist # allowed nt. mismatches between barcode and read
        # Check files exists before starting:
        self.AdapterRemoval_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['AdapterRemoval_dir'])
        for _, row in self.inp_file_df.iterrows(): # Pull out each merged fastq file
            basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
            merged_fastq_fn = '{}/{}.collapsed.bz2'.format(self.AdapterRemoval_dir_abs, basename)
            assert(os.path.exists(merged_fastq_fn))
        
        # Make output folder:
        self._make_dir(overwrite_dir)

    def _make_dir(self, overwrite):
        # Create folder for files:
        self.BC_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['BC_dir'])
        try:
            os.mkdir(self.BC_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.BC_dir_abs)
                os.mkdir(self.BC_dir_abs)
            else:
                print('Using existing folder because overwrite set to false: {}'.format(self.BC_dir_abs))

    def run_parallel(self, n_jobs=4, load_previous=False):
        '''
        Submit the merged files for barcode splitting.

        Keyword arguments:
        n_jobs -- Number of subprocesses to be started in parallel (default 4)
        load_previous -- Attempt to load results from a previous barcode split by looking up index-pair_stats.xlsx and sample_stats.xlsx (default False)
        '''
        # Must check that not too many file handles are opened at the same time:
        max_fh = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
        fh_div = len(self.sample_df)*3 // max_fh
        try:
            assert(fh_div == 0)
        except:
            raise Exception('Large number of samples, could cause problems with too many open filehandles under parallel processing. Either switch to serial processing, or implement parallel processing in smaller chunks.')

        if load_previous is False:
            # Run parallel:
            data = list(self.inp_file_df.iterrows())
            with WorkerPool(n_jobs=n_jobs) as pool:
                results = pool.map(self._split_file, data)
            self._collect_stats(results)
        else:
            try:
                self.inp_file_df = pd.read_excel('{}/index-pair_stats.xlsx'.format(self.BC_dir_abs), index_col=0)
                self.sample_df = pd.read_excel('{}/sample_stats.xlsx'.format(self.BC_dir_abs), index_col=0)
                print('Loaded results from previous run... Not running barcode split.')
            except Exception as err:
                print('Attempted to read previous stats from index-pair_stats and sample_stats, but failed...')
                raise err
        return(self.sample_df, self.inp_file_df)

    def _split_file(self, index, row):
        basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
        merged_fastq_fn = '{}/{}.collapsed.bz2'.format(self.AdapterRemoval_dir_abs, basename)
        # List the barcodes and associated sample names:
        mask = (self.sample_df['fastq_mate1_filename'] == row['fastq_mate1_filename']) & (self.sample_df['fastq_mate2_filename'] == row['fastq_mate2_filename'])
        assert(len({snu for snu in self.sample_df[mask]['sample_name_unique'].values}) == len(self.sample_df[mask]))
        assert(len({snu for snu in self.sample_df[mask]['barcode_seq'].values}) == len(self.sample_df[mask]))
        bc_fh = [(k, v, bz2.open('{}/{}.fastq.bz2'.format(self.BC_dir_abs, v), "wt")) for k, v in zip(self.sample_df[mask]['barcode_seq'].values, self.sample_df[mask]['sample_name_unique'].values)]
        unmapped_fh = bz2.open('{}/{}_unmapped.fastq.bz2'.format(self.BC_dir_abs, basename), "wt")

        # Collect stats:
        Nmapped = 0
        Nunmapped = 0
        Ncc = {k:0 for k in self.sample_df[mask]['sample_name_unique'].values}
        Ncca = {k:0 for k in self.sample_df[mask]['sample_name_unique'].values}
        Ntot = {k:0 for k in self.sample_df[mask]['sample_name_unique'].values}
        Nreadlen = {k: Counter() for k in self.sample_df[mask]['sample_name_unique'].values}

        # Iterate over each record in the fastq file:
        with bz2.open(merged_fastq_fn, "rt") as input_fh:
            for title, seq, qual in FastqGeneralIterator(input_fh):
                # Search for barcodes and write to barcode specific file:
                found = False
                for bc, sample_name, fh in bc_fh:
                    if ('N' in bc and sum(l1!=l2 for l1, l2 in zip(seq[-len(bc):], bc) if l2 != 'N') <= self.max_dist) or \
                       (not 'N' in bc and jellyfish.hamming_distance(seq[-len(bc):], bc) <= self.max_dist):
                    # Old if statement compatible with random nucleotides in barcode:
                    # if all(l1==l2 for l1, l2 in zip(seq[-len(bc):], bc) if l2 != 'N'):
                    # if jellyfish.hamming_distance(seq[-len(bc):], bc) <= self.max_dist:
                        found = True
                        # Add adapter sequence to title:
                        title = title + ':' + seq[-len(bc):]
                        fh.write("@{}\n{}\n+\n{}\n".format(title, seq[:-len(bc)], qual[:-len(bc)]))
                        Nmapped += 1
                        Ntot[sample_name] += 1
                        Nreadlen[sample_name][len(seq) - len(bc)] += 1
                        # Count if CC, CCA or not:
                        if seq[-(len(bc)+2):-len(bc)] == 'CC':
                            Ncc[sample_name] += 1
                        elif seq[-(len(bc)+3):-len(bc)] == 'CCA':
                            Ncca[sample_name] += 1
                        break
                if not found:
                    Nunmapped += 1
                    unmapped_fh.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
        for _, _, fh in bc_fh:
            fh.close()
        unmapped_fh.close()
        return(basename, Nmapped, Nunmapped, Ncc, Ncca, Ntot, Nreadlen)

    def _collect_stats(self, results):
        # Unfold the results output:
        Ncc_union = dict()
        Ncca_union = dict()
        Ntot_union = dict()
        Nmapped_dict = dict()
        Nunmapped_dict = dict()
        Nreadlen_union = dict()
        for res in results:
            basename, Nmapped, Nunmapped, Ncc, Ncca, Ntot, Nreadlen = res
            Nmapped_dict[basename] = Nmapped
            Nunmapped_dict[basename] = Nunmapped
            Ncc_union.update(Ncc)
            Ncca_union.update(Ncca)
            Ntot_union.update(Ntot)
            Nreadlen_union.update(Nreadlen)

        # Sort mapped/unmapped:
        Nmapped_list = list()
        Nunmapped_list = list()
        for _, row in self.inp_file_df.iterrows():
            basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
            Nmapped_list.append(Nmapped_dict[basename])
            Nunmapped_list.append(Nunmapped_dict[basename])

        # Add stats to input file info dataframe:
        self.inp_file_df['N_BC-mapped'] = Nmapped_list
        self.inp_file_df['N_BC-unmapped'] = Nunmapped_list
        self.inp_file_df['N_sum-check'] = self.inp_file_df['N_BC-mapped'] + self.inp_file_df['N_BC-unmapped']
        self.inp_file_df['percent_BC-mapped'] = self.inp_file_df['N_BC-mapped'].values / self.inp_file_df['N_merged'].values *100

        # Add stats to sample info dataframe:
        self.sample_df['N_total'] = [Ntot_union[sn] for sn in self.sample_df['sample_name_unique']]
        self.sample_df['N_CC'] = [Ncc_union[sn] for sn in self.sample_df['sample_name_unique']]
        self.sample_df['N_CCA'] = [Ncca_union[sn] for sn in self.sample_df['sample_name_unique']]
        self.sample_df['N_CCA+CC'] = self.sample_df['N_CCA'].values + self.sample_df['N_CC'].values
        # Allow division by zero, in case of small sample:
        with np.errstate(divide='ignore', invalid='ignore'):
            self.sample_df['CCA+CC_percent_total'] = self.sample_df['N_CCA+CC'].values / self.sample_df['N_total'].values *100
            self.sample_df['percent_CCA'] = self.sample_df['N_CCA'].values / self.sample_df['N_CCA+CC'].values *100
        # Sanity check:
        assert(self.inp_file_df['N_merged'].sum() - self.inp_file_df['N_BC-unmapped'].sum()) == self.sample_df['N_total'].sum()
        # Dump stats as Excel file:
        self.inp_file_df.to_excel('{}/index-pair_stats.xlsx'.format(self.BC_dir_abs))
        self.sample_df.to_excel('{}/sample_stats.xlsx'.format(self.BC_dir_abs))

        # Save read length distributions as CSV:
        readlen_rows = []
        for sn, counter in Nreadlen_union.items():
            for length, count in sorted(counter.items()):
                readlen_rows.append({'sample_name_unique': sn, 'read_length': length, 'count': count})
        if readlen_rows:
            readlen_df = pd.DataFrame(readlen_rows)
            readlen_df.to_csv('{}/read_length_distributions.csv'.format(self.BC_dir_abs), index=False)



class Kmer_analysis:
    '''
    Class to find Kmers at the end of unmapped reads,
    in order to determine if barcode mapping was efficient.

    Keyword arguments:
    k_size -- Kmer size (default 5)
    overwrite -- Overwrite previous Kmer analysis folder (default True)
    '''
    def __init__(self, dir_dict, inp_file_df, index_df, \
                 k_size=5, overwrite=True):
        # Input:
        self.inp_file_df, self.k_size = inp_file_df, k_size
        self.dir_dict = dir_dict

        self.index_dict = dict()
        for t, i, s in zip(index_df['type'].values, index_df['id'].values, \
                           index_df['sequence'].values):
            if t not in self.index_dict:
                self.index_dict[t] = dict()
            self.index_dict[t][i] = s

        self.filter_dict = dict()
        
        # Check files exists before starting:
        self.BC_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], \
                                            self.dir_dict['data_dir'], \
                                            self.dir_dict['BC_dir'])
        for _, row in self.inp_file_df.iterrows():
            basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
            unmapped_fn = '{}/{}_unmapped.fastq.bz2'.format(self.BC_dir_abs, basename)
            assert(os.path.exists(unmapped_fn))

        # Create folder for files:
        self.Kmer_dir_abs = '{}/{}'.format(self.BC_dir_abs, 'Kmer_analysis')
        try:
            os.mkdir(self.Kmer_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.Kmer_dir_abs)
                os.mkdir(self.Kmer_dir_abs)
            else:
                raise Exception('Folder exists and overwrite set to false: {}'.format(self.Kmer_dir_abs))

    def filter_3p_fasta(self, input_fasta, filter_search_size=7):
        '''
        Generate a filter composed of the Kmers contained
        in the last X nt. of a set of fasta sequences,
        e.g. all human tRNA seqeunces.
        '''
        with open(input_fasta, "r") as seq_fh:
            for seq_obj in SeqIO.parse(seq_fh, "fasta"):
                self.filter_dict = self._add_kmers(self.filter_dict, \
                                                    str(seq_obj.seq)[-filter_search_size:])

    def filter_window_BC(self, filter_window=(0, 11)):
        '''
        Generate a filter composed of the Kmers contained
        in a window of the input barcode sequences,
        e.g. constant region of an adapter.
        '''
        for bc, bc_seq in self.index_dict['barcode'].items():
            self.filter_dict = self._add_kmers(self.filter_dict, \
                                                bc_seq[filter_window[0]:filter_window[1]])

    def search_unmapped(self, search_size=13):
        '''Search for Kmers of size N=search_size'''
        # Find Kmers and store in dictionary:
        k_dict_all = dict()
        for _, row in self.inp_file_df.iterrows():
            k_dict = dict()
            basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
            unmapped_fn = '{}/{}_unmapped.fastq.bz2'.format(self.BC_dir_abs, basename)
            with bz2.open(unmapped_fn, "rt") as unmapped_fh:
                for title, seq, qual in FastqGeneralIterator(unmapped_fh):
                    if len(seq) >= search_size:
                        k_dict = self._add_kmers(k_dict, seq[-search_size:])
            
            self._write_stats(k_dict, basename)
            for kmer, obs in k_dict.items():
                try:
                    k_dict_all[kmer] += obs
                except KeyError:
                    k_dict_all[kmer] = obs

        kmer_df = self._write_stats(k_dict_all, 'ALL', return_df=True)
        return(kmer_df)

    def _write_stats(self, k_dict, outp_ext, return_df=False):
        # Rank Kmers by occurence and find closely related adapters: 
        kmer_df_dat = list()
        for kmer_seq, count in sorted(k_dict.items(), key=lambda x:x[1], reverse=True):
            bc_min_dist, dist_min = self._find_min_dist_bc(kmer_seq)
            if dist_min < 2:
                kmer_df_dat.append([kmer_seq, count, dist_min, bc_min_dist])
            else:
                kmer_df_dat.append([kmer_seq, count, None, None])
        kmer_df = pd.DataFrame(kmer_df_dat, columns=['Kmer', 'Count', \
                                                     'Barcode distance', 'Barcode'])
        kmer_df.to_excel('{}/{}_unmapped-Kmer-analysis.xlsx'.format(self.Kmer_dir_abs, outp_ext))
        if return_df:
            return(kmer_df)

    def _add_kmers(self, k_dict, seq):
        '''Find Kmers in input sequence and add to dictionary if not filtered.'''
        for i in range(len(seq) - self.k_size + 1):
            kmer = seq[i:(i+self.k_size)]
            if kmer not in self.filter_dict:
                try:
                    k_dict[kmer] += 1
                except KeyError:
                    k_dict[kmer] = 1
        return(k_dict)

    def _find_min_dist_bc(self, kmer_seq):
        '''Search for the Kmers in the adapter sequences.'''
        dist_min = 999
        bc_min_dist = ''
        for bc, bc_seq in self.index_dict['barcode'].items():
            for i in range(len(bc_seq) - len(kmer_seq) + 1):
                window = bc_seq[i:(i+len(kmer_seq))]
                dist = jellyfish.hamming_distance(window, kmer_seq)
                if dist < dist_min:
                    dist_min = dist
                    bc_min_dist = bc
        return(bc_min_dist, dist_min)



class BC_analysis:
    '''
    Class to find barcodes at the end of unmapped reads,
    in order to determine if barcode mapping was efficient.

    Keyword arguments:
    BC_size_3p -- Length of the barcode sequence to extract from the 3p end of each barcode in index_df (default 5)
    overwrite -- Overwrite previous barcode analysis folder (default True)
    '''
    def __init__(self, dir_dict, inp_file_df, index_df, \
                 BC_size_3p=5, overwrite=True):
        # Input:
        self.inp_file_df, self.BC_size_3p = inp_file_df, BC_size_3p
        self.dir_dict = dir_dict

        self.index_dict = dict()
        for t, i, s in zip(index_df['type'].values, index_df['id'].values, index_df['sequence'].values):
            if t not in self.index_dict:
                self.index_dict[t] = dict()
            self.index_dict[t][i] = s
        self.bc_list = [seq[-self.BC_size_3p:] for seq in self.index_dict['barcode'].values()]
        self.bc2name = {seq[-self.BC_size_3p:]: name for name, seq in self.index_dict['barcode'].items()}
        
        # Check files exists before starting:
        self.BC_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], \
                                            self.dir_dict['data_dir'], \
                                            self.dir_dict['BC_dir'])
        for _, row in self.inp_file_df.iterrows():
            basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
            unmapped_fn = '{}/{}_unmapped.fastq.bz2'.format(self.BC_dir_abs, basename)
            assert(os.path.exists(unmapped_fn))

        # Create folder for files:
        self.BCanalysis_dir_abs = '{}/{}'.format(self.BC_dir_abs, 'BC_analysis')
        try:
            os.mkdir(self.BCanalysis_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.BCanalysis_dir_abs)
                os.mkdir(self.BCanalysis_dir_abs)
            else:
                print('Using existing folder because overwrite set to false: {}'.format(self.BCanalysis_dir_abs))

    def search_unmapped(self, search_size=13, group_dist=1, \
                        load_previous=False):
        '''
        Submit search for barcodes in unmapped reads.
        Keyword arguments:
        search_size -- Distance from the 3p of the read to search within (default 13)
        group_dist -- Group distance smaller than or equal to in aggregated output file (default 1)
        load_previous -- Attempt to load results from a previous barcode search by looking up the X_unmapped-BC-analysis.xlsx results file (default False)
        '''
        if load_previous:
            self.sum_df = pd.read_excel('{}/{}{}_unmapped-BC-analysis.xlsx'.format(self.BCanalysis_dir_abs, 'ALL-groupby-dist-', group_dist), index_col=0)
            print('Loaded results from previous run... Not running barcode analysis.')
            return(self.sum_df)

        # Find closest barcode and store in dictionary:
        k_dict_all = {bc: {} for bc in self.bc_list}
        for _, row in self.inp_file_df.iterrows():
            k_dict = {bc: {} for bc in self.bc_list}
            basename = '-'.join(re.split('[_/]', row['fastq_mate1_filename'])[:-1])
            unmapped_fn = '{}/{}_unmapped.fastq.bz2'.format(self.BC_dir_abs, basename)
            with bz2.open(unmapped_fn, "rt") as unmapped_fh:
                for title, seq, qual in FastqGeneralIterator(unmapped_fh):
                    if len(seq) >= search_size:
                        seq_min_dist, bc_min_dist, dist_min = self._find_min_dist_bc(seq[-search_size:])
                        seq_dist = seq_min_dist + '-' + str(dist_min)
                        try:
                            k_dict[bc_min_dist][seq_dist] += 1
                        except KeyError:
                            k_dict[bc_min_dist][seq_dist] = 1
            # Write stats for file:
            self._write_stats(k_dict, basename)
            # Collect stats for all files:
            for bc in self.bc_list:
                for dist, obs in k_dict[bc].items():
                    try:
                        k_dict_all[bc][dist] += obs
                    except KeyError:
                        k_dict_all[bc][dist] = obs

        # Write stats for all files:
        self.bc_df = self._write_stats(k_dict_all, 'ALL', return_df=True)
        # Group by barcode name and make total count:
        mask = self.bc_df['Distance'] <= group_dist
        self.sum_df = self.bc_df.loc[mask, ['Name', 'Count']].groupby('Name').sum().sort_values(by=['Count'], ascending=False).reset_index()
        self.sum_df.to_excel('{}/{}{}_unmapped-BC-analysis.xlsx'.format(self.BCanalysis_dir_abs, 'ALL-groupby-dist-', group_dist))
        return(self.sum_df)

    def _write_stats(self, k_dict, outp_ext, return_df=False):
        # Rank barcodes by occurence and write stats to Excel file:
        kmer_df_dat = list()
        for bc in self.bc_list:
            for seq_dist, obs in k_dict[bc].items():
                target, dist = seq_dist.split('-')
                dist = int(dist)
                kmer_df_dat.append([self.bc2name[bc], bc, target, obs, dist])

        bc_df = pd.DataFrame(kmer_df_dat, columns=['Name', 'Barcode', 'Target', 'Count', 'Distance'])
        bc_df = bc_df.sort_values(by=['Distance', 'Count', 'Name'], ascending=[True, False, True]).reset_index(drop=True)
        bc_df.to_excel('{}/{}_unmapped-BC-analysis.xlsx'.format(self.BCanalysis_dir_abs, outp_ext))
        if return_df:
            return(bc_df)

    def _find_min_dist_bc(self, seq):
        '''Search for the closest barcode sequence.'''
        dist_min = 999
        seq_min_dist = ''
        bc_min_dist = ''
        for bc in self.bc_list:
            for i in range(len(seq) - self.BC_size_3p + 1):
                window = seq[i:(i+self.BC_size_3p)]
                dist = jellyfish.hamming_distance(window, bc)
                if dist < dist_min:
                    dist_min = dist
                    seq_min_dist = window
                    bc_min_dist = bc
        return(seq_min_dist, bc_min_dist, dist_min)



class UMI_trim:
    '''
    Class to trim off the UMI from each read,
    add it to the fasta header and generate statistics
    on the UMI use.
    After trimming it is possible to downsample the number
    of reads such that all samples have the same number of reads.

    Supports two trimming modes:
    - 'pyrimidine': Legacy mode. Takes first UMI_len nt, checks last nt is T/C.
    - 'anchored': Finds a constant anchor sequence after a variable-length stagger (0..max_stagger nt).
      The full prefix (stagger + anchor) is used as the UMI for deduplication.

    Keyword arguments:
    mode -- Trimming mode: 'anchored' or 'pyrimidine' (default 'anchored')
    anchor_seq -- Constant sequence to search for in anchored mode (required when mode='anchored')
    max_stagger -- Maximum random nt before the anchor sequence (default 3)
    anchor_max_dist -- Hamming distance tolerance for anchor matching (default 1)
    UMI_end -- Set of nucleotides possible on the last UMI position, pyrimidine mode only (default {'T', 'C'})
    overwrite_dir -- Overwrite previous UMI trim folder (default False)
    check_input -- Check if input files exist (default True)
    UMI_len -- UMI length, pyrimidine mode only (default 10)
    UMI_bins -- Number of possible UMIs (default depends on mode)
    downsample_absolute -- Maximum number of trimmed UMI sequences, randomly choose for downsampling (default 2e6)
    downsample_percentile -- Determine the 'downsample_absolute' variable as a percentile of the number of trimmed UMI sequences in all input samples (default False)
    '''
    def __init__(self, dir_dict, sample_df, mode='anchored',
                 anchor_seq=None, max_stagger=3, anchor_max_dist=1,
                 UMI_end={'T', 'C'},
                 overwrite_dir=False, check_input=True, UMI_len=10,
                 UMI_bins=None, downsample_absolute=2e6,
                 downsample_percentile=False):
        self.mode = mode
        if mode == 'anchored':
            if anchor_seq is None:
                raise ValueError("anchor_seq is required when mode='anchored'")
            valid_bases = set('ACGTacgt')
            if not all(c in valid_bases for c in anchor_seq):
                raise ValueError(f"anchor_seq contains invalid characters: {anchor_seq}")
            self.anchor_seq = anchor_seq.upper()
            self.max_stagger = max_stagger
            self.anchor_max_dist = anchor_max_dist
            # Number of possible UMIs: sum of 4^i for stagger lengths 0..max_stagger
            self.UMI_bins = UMI_bins if UMI_bins is not None else sum(4**i for i in range(max_stagger + 1))
        elif mode == 'pyrimidine':
            self.UMI_len = UMI_len
            self.UMI_end = UMI_end
            # 9x random nt (A/G/T/C) and one purine (A/G) → 4^9 * 2
            self.UMI_bins = UMI_bins if UMI_bins is not None else 4**9 * 2
        else:
            raise ValueError(f"Unknown UMI trim mode: {mode!r}. Use 'anchored' or 'pyrimidine'.")
        
        # Input:
        self.sample_df = sample_df
        self.dir_dict = dir_dict
        self.downsample_absolute = downsample_absolute
        self.downsample_percentile = downsample_percentile

        # Check files exists before starting:
        self.BC_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['BC_dir'])
        if check_input:
            for _, row in self.sample_df.iterrows():
                mapped_fn = '{}/{}.fastq.bz2'.format(self.BC_dir_abs, row['sample_name_unique'])
                assert(os.path.exists(mapped_fn))

        # Make output folder:
        self._make_dir(overwrite_dir)

    def _make_dir(self, overwrite):
        # Create folder for files:
        self.UMI_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['UMI_dir'])
        try:
            os.mkdir(self.UMI_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.UMI_dir_abs)
                os.mkdir(self.UMI_dir_abs)
            else:
                print('Using existing folder because overwrite set to false: {}'.format(self.UMI_dir_abs))

    def run_parallel(self, n_jobs=4, load_previous=False):
        '''
        Submit files for UMI trimming.

        Keyword arguments:
        n_jobs -- Number of subprocesses to be started in parallel (default 4)
        load_previous -- Attempt to load results from a previous UMI trim by looking up sample_stats.xlsx (default False)
        '''
        if load_previous is False:
            # Run parallel:
            data = list(self.sample_df.iterrows())
            with WorkerPool(n_jobs=n_jobs) as pool:
                results = pool.map(self._trim_file, data)

            # Find the max number of reads for downsampling:
            ds_abs_max = None
            if self.downsample_percentile:
                stats_df = pd.DataFrame(results, columns=['sample_name_unique', 'N_after_trim', \
                                                          'N_UMI_observed', 'N_UMI_expected', \
                                                          'N_after_downsample'])
                ds_abs_max = int(np.percentile(stats_df['N_after_trim'], self.downsample_percentile))
            elif self.downsample_absolute:
                ds_abs_max = int(self.downsample_absolute)
            # Perform downsampling if requested:
            if not ds_abs_max is None:
                print('Downsampling UMI trimmed sequences to maximum {} reads.'.format(ds_abs_max))
                # Update results with the downsampling choice:
                for ri in range(len(results)):
                    if results[ri][-1] > ds_abs_max:
                        results[ri][-1] = ds_abs_max
                # Run downsampling in parallel:
                with WorkerPool(n_jobs=n_jobs) as pool:
                    _ = pool.map(self._downsample_file, results)

            self._collect_stats(results)
        else:
            try:
                self.sample_df = pd.read_excel('{}/sample_stats.xlsx'.format(self.UMI_dir_abs), index_col=0)
                print('Loaded results from previous run... Not running UMI trimming.')
            except Exception as err:
                print('Attempted to read previous stats from sample_stats, but failed...')
                raise err
        return(self.sample_df)

    def _trim_file(self, index, row):
        input_fnam = '{}/{}.fastq.bz2'.format(self.BC_dir_abs, row['sample_name_unique'])
        output_fnam = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
        output_fnam_untrimmed = '{}/{}_untrimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
        untrimmed_fh = bz2.open(output_fnam_untrimmed, "wt")

        UMIs = set()
        Nseqs = 0
        with bz2.open(output_fnam, "wt") as output_fh:
            with bz2.open(input_fnam, "rt") as input_fh:
                if self.mode == 'pyrimidine':
                    for title, seq, qual in FastqGeneralIterator(input_fh):
                        umi = seq[0:self.UMI_len]
                        if umi[-1] in self.UMI_end: # UMI found
                            UMIs.add(umi)
                            Nseqs += 1
                            title = title + ':' + umi
                            output_fh.write("@{}\n{}\n+\n{}\n".format(title, seq[self.UMI_len:], qual[self.UMI_len:]))
                        else: # UMI not found
                            untrimmed_fh.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
                elif self.mode == 'anchored':
                    anchor_len = len(self.anchor_seq)
                    for title, seq, qual in FastqGeneralIterator(input_fh):
                        found = False
                        for stagger in range(self.max_stagger + 1):
                            if stagger + anchor_len > len(seq):
                                break
                            candidate = seq[stagger:stagger + anchor_len]
                            if jellyfish.hamming_distance(candidate, self.anchor_seq) <= self.anchor_max_dist:
                                trim_pos = stagger + anchor_len
                                umi = seq[0:trim_pos]
                                UMIs.add(umi)
                                Nseqs += 1
                                title = title + ':' + umi
                                output_fh.write("@{}\n{}\n+\n{}\n".format(title, seq[trim_pos:], qual[trim_pos:]))
                                found = True
                                break
                        if not found:
                            untrimmed_fh.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
        untrimmed_fh.close()

        # Calculate and return the observed and expected UMI count:
        N_umi_obs = len(UMIs)
        N_umi_exp = self.UMI_bins*(1-((self.UMI_bins-1) / self.UMI_bins)**Nseqs)
        return([row['sample_name_unique'], Nseqs, N_umi_obs, N_umi_exp, Nseqs])

    def _downsample_file(self, unam, Nseqs, unused1, unused2, Nds):
        if Nseqs == Nds:
            return()

        # Randomly sample indices:
        idx_set = set(np.random.choice(Nseqs, size=Nds, replace=False))
        # I/O and read sampling:
        UMI_fnam = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, unam)
        UMI_fnam_DS = '{}/{}_UMI-trimmed.fastq.bz2_DStmp'.format(self.UMI_dir_abs, unam)
        with bz2.open(UMI_fnam, "rt") as fh_in:
            with bz2.open(UMI_fnam_DS, "wt") as fh_out:
                for idx, fq_in in enumerate(FastqGeneralIterator(fh_in)):
                    if idx in idx_set:
                        fh_out.write("@{}\n{}\n+\n{}\n".format(*fq_in))
        # Overwrite old UMI file with new:
        os.replace(UMI_fnam_DS, UMI_fnam)

    def _collect_stats(self, results):
        # Stats to dataframe:
        stats_df = pd.DataFrame(results, columns=['sample_name_unique', 'N_after_trim', \
                                                  'N_UMI_observed', 'N_UMI_expected', \
                                                  'N_after_downsample'])
        # Merge stats with sample info dataframe:
        self.sample_df = self.sample_df.drop(columns=['N_after_trim', 'N_UMI_observed', 'N_UMI_expected', 'N_after_downsample'], errors='ignore')
        self.sample_df = self.sample_df.merge(stats_df, on=['sample_name_unique'])        
        # Add stats:
        self.sample_df['percent_seqs_after_UMI_trim'] = self.sample_df['N_after_trim'] / self.sample_df['N_total'] * 100
        self.sample_df['percent_UMI_obs-vs-exp'] = self.sample_df['N_UMI_observed'] / self.sample_df['N_UMI_expected'] * 100
        # Rearrange column order:
        colnames = self.sample_df.columns.tolist()
        new_colnames = colnames[:-3] + colnames[-2:] + [colnames[-3]]
        self.sample_df = self.sample_df.loc[:, new_colnames].copy()
        # Dump stats as Excel file:
        self.sample_df.to_excel('{}/sample_stats.xlsx'.format(self.UMI_dir_abs))



class _5p_trim:
    '''
    Class to trim off 5-prime adapter from each read.

    Keyword arguments:
    overwrite_dir -- Overwrite previous UMI trim folder (default False)
    check_input -- Check if input files exist (default True)
    downsample_absolute -- Maximum number of trimmed UMI sequences, randomly choose for downsampling (default 2e6)
    downsample_percentile -- Determine the 'downsample_absolute' variable as a percentile of the number of trimmed UMI sequences in all input samples (default False)
    '''
    def __init__(self, dir_dict, sample_df, adapters, \
                 overwrite_dir=False, check_input=True, \
                 downsample_absolute=2e6, \
                 downsample_percentile=False, \
                 max_dist=0):
        self.adapters = adapters
        self.max_dist = max_dist

        # Input:
        self.sample_df = sample_df
        self.dir_dict = dir_dict
        self.downsample_absolute = downsample_absolute
        self.downsample_percentile = downsample_percentile

        # Check files exists before starting:
        self.BC_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['BC_dir'])
        if check_input:
            for _, row in self.sample_df.iterrows():
                mapped_fn = '{}/{}.fastq.bz2'.format(self.BC_dir_abs, row['sample_name_unique'])
                assert(os.path.exists(mapped_fn))

        # Make output folder:
        self._make_dir(overwrite_dir)

    def _make_dir(self, overwrite):
        # Create folder for files:
        self.UMI_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['UMI_dir'])
        try:
            os.mkdir(self.UMI_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.UMI_dir_abs)
                os.mkdir(self.UMI_dir_abs)
            else:
                print('Using existing folder because overwrite set to false: {}'.format(self.UMI_dir_abs))

    def run_parallel(self, n_jobs=4, load_previous=False):
        '''
        Submit files for trimming.

        Keyword arguments:
        n_jobs -- Number of subprocesses to be started in parallel (default 4)
        load_previous -- Attempt to load results from a previous trim by looking up sample_stats.xlsx (default False)
        '''
        if load_previous is False:
            # Run parallel:
            data = list(self.sample_df.iterrows())
            with WorkerPool(n_jobs=n_jobs) as pool:
                results = pool.map(self._trim_file, data)

            # Find the max number of reads for downsampling:
            ds_abs_max = None
            if self.downsample_percentile:
                stats_df = pd.DataFrame(results, columns=['sample_name_unique', 'N_after_trim', \
                                                          'N_after_downsample'])
                ds_abs_max = int(np.percentile(stats_df['N_after_trim'], self.downsample_percentile))
            elif self.downsample_absolute:
                ds_abs_max = int(self.downsample_absolute)
            # Perform downsampling if requested:
            if not ds_abs_max is None:
                print('Downsampling trimmed sequences to maximum {} reads.'.format(ds_abs_max))
                # Update results with the downsampling choice:
                for ri in range(len(results)):
                    if results[ri][-1] > ds_abs_max:
                        results[ri][-1] = ds_abs_max
                # Run downsampling in parallel:
                with WorkerPool(n_jobs=n_jobs) as pool:
                    _ = pool.map(self._downsample_file, results)

            self._collect_stats(results)
        else:
            try:
                self.sample_df = pd.read_excel('{}/sample_stats.xlsx'.format(self.UMI_dir_abs), index_col=0)
                print('Loaded results from previous run... Not running trimming.')
            except Exception as err:
                print('Attempted to read previous stats from sample_stats, but failed...')
                raise err
        return(self.sample_df)

    def _trim_file(self, index, row):
        input_fnam = '{}/{}.fastq.bz2'.format(self.BC_dir_abs, row['sample_name_unique'])
        output_fnam = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
        output_fnam_untrimmed = '{}/{}_untrimmed.fastq.bz2'.format(self.UMI_dir_abs, row['sample_name_unique'])
        untrimmed_fh = bz2.open(output_fnam_untrimmed, "wt")

        Nseqs = 0
        with bz2.open(output_fnam, "wt") as output_fh:
            with bz2.open(input_fnam, "rt") as input_fh:
                for title, seq, qual in FastqGeneralIterator(input_fh):
                    # Search for adapters:
                    found = False
                    for bc in self.adapters:
                        if ('N' in bc and sum(l1!=l2 for l1, l2 in zip(seq[:len(bc)], bc) if l2 != 'N') <= self.max_dist) or \
                           (not 'N' in bc and jellyfish.hamming_distance(seq[:len(bc)], bc) <= self.max_dist):
                            found = True
                            # Add adapter sequence to title:
                            title = title + ':' + seq[:len(bc)]
                            output_fh.write("@{}\n{}\n+\n{}\n".format(title, seq[len(bc):], qual[len(bc):]))
                            Nseqs += 1
                            break
                    if not found:
                        untrimmed_fh.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
        untrimmed_fh.close()
        return([row['sample_name_unique'], Nseqs, Nseqs])

    def _downsample_file(self, unam, Nseqs, Nds):
        if Nseqs == Nds:
            return()

        # Randomly sample indices:
        idx_set = set(np.random.choice(Nseqs, size=Nds, replace=False))
        # I/O and read sampling:
        UMI_fnam = '{}/{}_UMI-trimmed.fastq.bz2'.format(self.UMI_dir_abs, unam)
        UMI_fnam_DS = '{}/{}_UMI-trimmed.fastq.bz2_DStmp'.format(self.UMI_dir_abs, unam)
        with bz2.open(UMI_fnam, "rt") as fh_in:
            with bz2.open(UMI_fnam_DS, "wt") as fh_out:
                for idx, fq_in in enumerate(FastqGeneralIterator(fh_in)):
                    if idx in idx_set:
                        fh_out.write("@{}\n{}\n+\n{}\n".format(*fq_in))
        # Overwrite old UMI file with new:
        os.replace(UMI_fnam_DS, UMI_fnam)

    def _collect_stats(self, results):
        # Stats to dataframe:
        stats_df = pd.DataFrame(results, columns=['sample_name_unique', 'N_after_trim', \
                                                  'N_after_downsample'])
        # Merge stats with sample info dataframe:
        self.sample_df = self.sample_df.drop(columns=['N_after_trim', 'N_after_downsample'], errors='ignore')
        self.sample_df = self.sample_df.merge(stats_df, on=['sample_name_unique'])
        # Add stats:
        self.sample_df['percent_seqs_after_UMI_trim'] = self.sample_df['N_after_trim'] / self.sample_df['N_total'] * 100
        # Rearrange column order:
        # colnames = self.sample_df.columns.tolist()
        # new_colnames = colnames[:-3] + colnames[-2:] + [colnames[-3]]
        # self.sample_df = self.sample_df.loc[:, new_colnames].copy()
        # Dump stats as Excel file:
        self.sample_df.to_excel('{}/sample_stats.xlsx'.format(self.UMI_dir_abs))



