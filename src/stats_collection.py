import os, shutil, bz2, warnings, json, gc
import json_stream
import pysam
from Bio import SeqIO
import pandas as pd
import numpy as np
import csv  # Add to the imports at the top if it's not there already

from mpire import WorkerPool
import warnings



class STATS_collection:
    '''
    Class to collect statistics from the
    alignment results.

    Keyword arguments:
    common_seqs -- bzip2 compressed fasta file of commonly observed sequences to avoid duplicated alignments (default None)
    ignore_common_count -- Ignore common count even if X_common-seq-obs.json filename exists (default False)
    check_exists -- Check if required files exist before starting (default True)
    overwrite_dir -- Overwrite old stats folder if any exists (default False)
    reads_SW_sorted -- Assume reads and SW results are sorted in the same order. This gives a massive memory saving. (default True)
    from_UMIdir -- Is the input data from a folder made with the UMI_trim class? (default True)
    max_5p_non_temp -- Maximum length of 5p non templated nucleotides (default 10)
    UMI_bins -- Number of possible UMIs (default 4^9 x 2)
    '''
    def __init__(self, dir_dict, tRNA_data, sample_df, common_seqs=None, \
                 ignore_common_count=False, check_exists=True, overwrite_dir=False, \
                 reads_SW_sorted=True, from_UMIdir=True, max_5p_non_temp=10, \
                 UMI_bins=4**9*2):

        # Update the CSV Headers
        self.stats_csv_header = ['readID', 'common_seq', 'sample_name_unique', \
                                 'sample_name', 'replicate', 'barcode', 'species', 'tRNA_annotation', \
                                 'align_score', 'fmax_score', 'Ndeletions', 'Ninsertions', \
                                 'unique_annotation', 'tRNA_annotation_len', \
                                 'align_5p_idx', 'align_3p_idx', 'align_5p_nt', 'align_3p_nt', \
                                 'codon', 'anticodon', 'amino_acid', '5p_cover', '3p_cover', \
                                 '5p_non-temp', '3p_non-temp', '5p_UMI', '3p_BC', \
                                 'align_gap', 'fmax_score>0.9', 'UMIcount', 'count']

        # Adding new classification fields
        self.stats_csv_header.extend([
            "is_5p_tRF", "is_degraded_tRNA", "is_pre_tRNA",
            "5p_non_temp_seq", "3p_non_temp_seq"
        ])

        # Update Aggregated Statistics Columns
        self.stats_agg_cols = ['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'species', \
                               'tRNA_annotation', 'tRNA_annotation_len', 'unique_annotation', \
                               '5p_cover', '3p_cover', 'align_3p_nt', 'codon', 'anticodon', 'amino_acid', \
                               'align_gap', 'fmax_score>0.9', 'UMIcount', 'count', 'is_5p_tRF', \
                               'is_degraded_tRNA', 'is_pre_tRNA']

        self.tRNA_data, self.sample_df = tRNA_data, sample_df
        self.dir_dict = dir_dict
        self.common_seqs_fnam = common_seqs
        self.reads_SW_sorted = reads_SW_sorted
        self.from_UMIdir = from_UMIdir
        self.max_5p_non_temp = max_5p_non_temp
        self.UMI_bins = UMI_bins

        # Check files exist before starting:
        self.align_dir_abs = f"{self.dir_dict['NBdir']}/{self.dir_dict['data_dir']}/{self.dir_dict['align_dir']}"
        if self.from_UMIdir:
            self.UMI_dir_abs = f"{self.dir_dict['NBdir']}/{self.dir_dict['data_dir']}/{self.dir_dict['UMI_dir']}"

        if check_exists:
            for _, row in self.sample_df.iterrows():
                SWres_fnam = f"{self.align_dir_abs}/{row['sample_name_unique']}_SWalign.json.bz2"
                assert os.path.exists(SWres_fnam)

                if self.from_UMIdir:
                    trimmed_fn = f"{self.UMI_dir_abs}/{row['sample_name_unique']}_UMI-trimmed.fastq.bz2"
                    assert os.path.exists(trimmed_fn)
                else:
                    fpath = f"{self.dir_dict['NBdir']}/{self.dir_dict['data_dir']}/{row['path']}"
                    assert os.path.exists(fpath)

                common_obs_fn = f"{self.align_dir_abs}/{row['sample_name_unique']}_common-seq-obs.json"
                if self.common_seqs_fnam:
                    assert os.path.exists(common_obs_fn)
                elif not self.common_seqs_fnam and not ignore_common_count and os.path.exists(common_obs_fn):
                    raise Exception(f"Found common sequence counts for {row['sample_name_unique']}, but common sequences are not specified.")

        # Create output directory
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

    def run_parallel(self, n_jobs=4, verbose=True, load_previous=False):
        if load_previous:
            try:
                stats_agg_fnam = f"{self.stats_dir_abs}/ALL_stats_aggregate.csv"
                self.concat_df = pd.read_csv(stats_agg_fnam, keep_default_na=False, dtype=self.stats_agg_cols_td)
                print("Loaded previous results.")
                return self.concat_df
            except Exception as err:
                print(f"Error loading previous results: {err}")

        print("Starting parallel processing...")
        
        data = list(self.sample_df.iterrows())

        # Use a controlled chunk size to prevent memory overflow
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self._collect_stats, data, chunk_size=10)

        self._concat_stats(results)
        return self.concat_df

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
        try:
            print(f"Processing sample: {row['sample_name_unique']}")

            stats_fnam = f"{self.stats_dir_abs}/{row['sample_name_unique']}_stats.csv.bz2"
            stats_agg_fnam = f"{self.stats_dir_abs}/{row['sample_name_unique']}_stats_aggregate.csv"

            # Ensure file exists before proceeding
            if not os.path.exists(stats_fnam):
                raise FileNotFoundError(f"Missing required file: {stats_fnam}")

            # Read compressed CSV safely
            with bz2.open(stats_fnam, 'rt') as stats_fh:
                stat_df = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)

            # Define 5' tRNA Fragments, Degraded tRNAs, and Pre-tRNAs
            stat_df["is_5p_tRF"] = (stat_df["align_5p_idx"] == 1) & \
                                (stat_df["align_3p_idx"] < stat_df["tRNA_annotation_len"]) & \
                                (stat_df["3p_non-temp"].str.len() == 0)

            stat_df["is_degraded_tRNA"] = (stat_df["align_5p_idx"] > 1) & \
                                        (stat_df["align_3p_idx"] < stat_df["tRNA_annotation_len"])

            # Pre-tRNA (combining 5' and 3' unprocessed)
            stat_df["is_pre_tRNA"] = (stat_df["5p_non-temp"].str.len() > 0) | (stat_df["3p_non-temp"].str.len() > 0)

            # Store 5' and 3' non-templated sequences for RT artifact checks
            stat_df["5p_non_temp_seq"] = stat_df["5p_non-temp"]
            stat_df["3p_non_temp_seq"] = stat_df["3p_non-temp"]

            # Filtering: Only Keep Full-Length tRNAs for Charge State Analysis
            row_mask = stat_df["align_3p_idx"] == stat_df["tRNA_annotation_len"]  # Full-length tRNAs only

            filtered_stat_df = stat_df[row_mask].copy()

            # Aggregate Information: Classify Reads Based on 3' End (CA, GA, CC, CG Separately)
            agg_df = filtered_stat_df.groupby(self.stats_agg_cols[:-2], as_index=False).agg(
                {"count": "sum", "UMIcount": "sum"}
            )

            # Separate Each 3' End Category (CA, GA, CC, CG)
            CA_mask = filtered_stat_df["align_3p_nts"] == "CA"
            GA_mask = filtered_stat_df["align_3p_nts"] == "GA"
            CC_mask = filtered_stat_df["align_3p_nts"] == "CC"
            CG_mask = filtered_stat_df["align_3p_nts"] == "CG"

            CA_agg = filtered_stat_df[CA_mask].groupby("tRNA_annotation").agg(Total_CA=("count", "sum")).reset_index()
            GA_agg = filtered_stat_df[GA_mask].groupby("tRNA_annotation").agg(Total_GA=("count", "sum")).reset_index()
            CC_agg = filtered_stat_df[CC_mask].groupby("tRNA_annotation").agg(Total_CC=("count", "sum")).reset_index()
            CG_agg = filtered_stat_df[CG_mask].groupby("tRNA_annotation").agg(Total_CG=("count", "sum")).reset_index()

            # Merge with aggregated tRNA counts
            agg_df = agg_df.merge(CA_agg, on="tRNA_annotation", how="left").fillna(0)
            agg_df = agg_df.merge(GA_agg, on="tRNA_annotation", how="left").fillna(0)
            agg_df = agg_df.merge(CC_agg, on="tRNA_annotation", how="left").fillna(0)
            agg_df = agg_df.merge(CG_agg, on="tRNA_annotation", how="left").fillna(0)

            # Add 5' Fragments, Degraded tRNA, and Pre-tRNA Categories
            tRNA_fragments = stat_df.groupby("tRNA_annotation").agg(
                Total_5p_tRFs=("is_5p_tRF", "sum"),
                Total_Degraded=("is_degraded_tRNA", "sum"),
                Total_Pre_tRNA=("is_pre_tRNA", "sum"),
            ).reset_index()

            agg_df = agg_df.merge(tRNA_fragments, on="tRNA_annotation", how="left").fillna(0)

            # Save Final Aggregated Results
            agg_df.to_csv(stats_agg_fnam, index=False)
            print(f"Saved aggregated stats to {stats_agg_fnam}")

            return agg_df

        except Exception as e:
            print(f"Error in _collect_stats for sample {row['sample_name_unique']}: {e}")
            return None

    def _read_non_common(self, row, stats_fh):
        import logging
        logging.basicConfig(filename='{}/{}/error_read_non_common.log'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir']), level=logging.ERROR)

        # Prepare the debug CSV file path to log tRNA_annotation and tRNA_annotation_len
        debug_csv_path = '{}/{}/tRNA_annotation_debug.csv'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'])

        # Open the debug CSV file for writing outside the loop to keep it open throughout the process
        debug_csv_file = open(debug_csv_path, 'a', newline='')  # Open in append mode
        csv_writer = csv.writer(debug_csv_file)

        # Write header only if the file is new or empty
        if debug_csv_file.tell() == 0:
            csv_writer.writerow(["tRNA_annotation", "tRNA_annotation_len", "align_3p_idx", "align_3p_nts"])

        try:
            # Read fastq files must be read to extract UMI and 5/3p non-template bases:
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
                # Parse JSON data as a stream (saves memory)
                if self.stream:
                    SWres = json_stream.load(SWres_fh)
                    SWres_json = SWres.persistent()
                else:
                    SWres_json = json.load(SWres_fh)

                # Loop through each read in the alignment results:
                for SWreadID, align_dict in SWres_json.items():
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
                    except Exception as e:
                        logging.error(f'Read ID ({SWreadID}) not found among read sequences in ({trimmed_fn}). Error: {e}')
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

                    # Get the tRNA annotation length
                    tRNA_annotation_len = self.tRNA_data[tRNA_annotation_first]['len']

                    align_5p_idx, align_3p_idx = align_dict['dpos']
                    align_5p_nt = align_dict['qseq'][0]
                    align_3p_nts = align_dict['qseq'][-2:]
                    
                    # Log tRNA_annotation and tRNA_annotation_len to CSV
                    csv_writer.writerow([tRNA_annotation, tRNA_annotation_len, align_3p_idx, align_3p_nts])
                    
                    qpos = align_dict['qpos']
                    _5p_non_temp = read_seq[0:(qpos[0]-1)]
                    _3p_non_temp = read_seq[qpos[1]:]
                    
                    # This is a special case for reads with 3' cleaved A for a 3'CCA tailed tRNA (because the SWAlign with a masked database doesn't align these reads right up to the mask) 
                    if align_3p_idx == (tRNA_annotation_len - 2) and _3p_non_temp == 'C':
                        align_3p_idx += 2
                        align_3p_nts = 'CC'
                        _3p_non_temp = ''

                    # This is a special case for reads with 3' cleaved A for a 3'CGA tailed tRNA (because the SWAlign with a masked database doesn't align these reads right up to the mask) 
                    if align_3p_idx == (tRNA_annotation_len - 2) and _3p_non_temp == 'G':
                        align_3p_idx += 2
                        align_3p_nts = 'CG'
                        _3p_non_temp = ''

                    # Move index for reads with 3' cleaved A:
                    if align_3p_idx == (tRNA_annotation_len - 1) and (align_3p_nts == 'CC' or align_3p_nts == 'CG'):
                        align_3p_idx += 1
                    
                    _3p_bc = row_exist_or_none(row, 'barcode_seq')
                    codon = self.tRNA_data[tRNA_annotation_first]['codon']
                    anticodon = self.tRNA_data[tRNA_annotation_first]['anticodon']
                    amino_acid = self.tRNA_data[tRNA_annotation_first]['amino_acid']
                    _5p_cover = align_5p_idx == 1
                    _3p_cover = align_3p_idx == tRNA_annotation_len

                    # Make booleans for gap and align score:
                    align_gap = Ndel > 0 or Nins > 0
                    fmax_score_09 = fmax_score > 0.9

                    # For "non-common" sequences multiple reads have not been collapsed:
                    UMIcount = 1
                    count = 1

                    # Print line to output csv file:
                    line_lst = [readID, common_seq, sample_name_unique, sample_name, replicate, \
                                barcode, species, tRNA_annotation, align_score, fmax_score, \
                                Ndel, Nins, unique_annotation, \
                                tRNA_annotation_len, align_5p_idx, align_3p_idx, align_5p_nt, \
                                align_3p_nts, codon, anticodon, amino_acid, _5p_cover, _3p_cover, \
                                _5p_non_temp, _3p_non_temp, _5p_umi, _3p_bc, \
                                align_gap, fmax_score_09, UMIcount, count]
                    csv_line = ','.join(map(str, line_lst))
                    print(csv_line, file=stats_fh)

            reads_fh.close()
        except Exception as e:
            logging.error(f"Error during CSV writing: {e}")
    
        finally:
            # Ensure the debug CSV file is closed after processing
            debug_csv_file.close()

    def _read_common(self, row, stats_fh):
        import logging
        logging.basicConfig(filename='{}/{}/error_read_common.log'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir']), level=logging.ERROR)

        # Prepare the debug CSV file path to log tRNA_annotation and tRNA_annotation_len
        debug_csv_path = '{}/{}/tRNA_annotation_debug.csv'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'])

        # Open the debug CSV file for writing outside the loop to keep it open
        debug_csv_file = open(debug_csv_path, 'a', newline='')  # Open in append mode
        csv_writer = csv.writer(debug_csv_file)

        # Write header only if the file is new or empty
        if debug_csv_file.tell() == 0:
            csv_writer.writerow(["tRNA_annotation", "tRNA_annotation_len", "align_3p_idx", "align_3p_nts"])
        
        try:
            # Read common sequences observations for this sample:
            common_obs_fn = '{}/{}_common-seq-obs.json'.format(self.align_dir_abs, row['sample_name_unique'])
            with open(common_obs_fn, 'r') as fh_in:
                obs_UMI_json = json.load(fh_in)
                common_obs = obs_UMI_json['common_obs']
                UMI_obs = obs_UMI_json['UMI_obs']

            # Open the alignment results:
            SWres_fnam = '{}/{}_SWalign.json.bz2'.format(self.align_dir_abs, 'common-seqs')
            with bz2.open(SWres_fnam, 'rt', encoding="utf-8") as SWres_fh:
                # Parse JSON data as a stream (saves memory)
                if self.stream:
                    SWres = json_stream.load(SWres_fh)
                    SWres_json = SWres.persistent()
                else:
                    SWres_json = json.load(SWres_fh)

                # Loop through each read in the alignment results:
                for readID, align_dict in SWres_json.items():
                    common_seq = True
                    readID_int = int(readID)
                    
                    # Skip reads that were not aligned:
                    if not align_dict['aligned']:
                        continue
                    # Skip if this sample did not have this common sequence:
                    elif common_obs[readID_int] == 0:
                        continue

                    # Collect the information:
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

                    # Get the tRNA annotation length
                    tRNA_annotation_len = self.tRNA_data[tRNA_annotation_first]['len']



                    align_5p_idx, align_3p_idx = align_dict['dpos']
                    align_5p_nt = align_dict['qseq'][0]
                    align_3p_nts = align_dict['qseq'][-2:]
                    # Log tRNA_annotation and tRNA_annotation_len to the CSV
                    csv_writer.writerow([tRNA_annotation, tRNA_annotation_len, align_3p_idx, align_3p_nts])

                    # Move index for reads with 3' cleaved A:
                    if align_3p_idx == (tRNA_annotation_len - 1) and align_3p_nts == 'CC':
                        align_3p_idx += 1

                    codon = self.tRNA_data[tRNA_annotation_first]['codon']
                    anticodon = self.tRNA_data[tRNA_annotation_first]['anticodon']
                    amino_acid = self.tRNA_data[tRNA_annotation_first]['amino_acid']
                    _5p_cover = align_5p_idx == 1
                    _3p_cover = align_3p_idx == tRNA_annotation_len

                    # Extract non-template bases from common reads:
                    try:
                        seq = self.common_seqs_info[readID]
                    except KeyError as e:
                        logging.error(f'Read ID ({readID}) not found in ({SWres_fnam}) among sequences. Error: {e}')
                        continue
                    except Exception as e:
                        logging.error(f'Error processing read {readID}: {e}')
                        continue

                    qpos = align_dict['qpos']
                    _5p_non_temp = seq[0:(qpos[0]-1)]
                    _3p_non_temp = seq[qpos[1]:]
                    _5p_umi = ''  # UMI information is lost when using common sequences
                    _3p_bc = row_exist_or_none(row, 'barcode_seq')

                    # Make booleans for gap and align score:
                    align_gap = Ndel > 0 or Nins > 0
                    fmax_score_09 = fmax_score > 0.9

                    # For common sequences add the read count:
                    UMIcount = int(UMI_obs[readID_int])
                    count = int(common_obs[readID_int])

                    # Print line to output csv file:
                    line_lst = [readID, common_seq, sample_name_unique, sample_name, replicate, \
                                barcode, species, tRNA_annotation, align_score, fmax_score, \
                                Ndel, Nins, unique_annotation, \
                                tRNA_annotation_len, align_5p_idx, align_3p_idx, align_5p_nt, \
                                align_3p_nts, codon, anticodon, amino_acid, _5p_cover, _3p_cover, \
                                _5p_non_temp, _3p_non_temp, _5p_umi, _3p_bc, \
                                align_gap, fmax_score_09, UMIcount, count]
                    csv_line = ','.join(map(str, line_lst))
                    print(csv_line, file=stats_fh)
        
        except Exception as e:
            logging.error(f"Error during CSV writing: {e}")

        finally:
            # Ensure the debug CSV file is closed after processing
            debug_csv_file.close()

    def _concat_stats(self, csv_paths):
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
    
    def make_bam_from_json(self, tRNA_list, bam_prefix='', overwrite=False):
        """
        Parses a JSON file and creates a BAM file for a specified tRNA annotation.

        Parameters:
        - json_file (str): Path to the JSON file containing alignment data.
        - tRNA_annotation (str): The tRNA annotation to filter and create the BAM file for.
        - output_dir (str): Directory to save the BAM file. Defaults to stats directory.
        - bam_prefix (str): Prefix for the BAM file name. Defaults to "tRNA".

        Returns:
        - str: Path to the created BAM file.
        """
        """
        Create BAM files for all samples in sample_df and a specified tRNA_annotation.

        Args:
            tRNA_annotation (str): tRNA annotation to filter reads for the BAM file.
            bam_prefix (str): Optional prefix for BAM filenames. Default is ''.
            overwrite (bool): Whether to overwrite existing BAM files. Default is False.
        """
        from pysam import AlignmentFile
        # Let's create a bam file directory
        bam_dir = f"{self.align_dir_abs}/bam_files"

        # Iterate over sample_df to extract sample_name_unique
        for _, row in self.sample_df.iterrows():
            sample_name = row['sample_name_unique']

            # Define the input JSON file path
            json_file = f"{self.align_dir_abs}/{sample_name}_SWalign.json.bz2"

            # Loop over tRNA_list to create BAM files for each tRNA_annotation
            for tRNA_annotation in tRNA_list:
                # Define the output BAM file path
                bam_filename = f"{bam_prefix}{'_' if bam_prefix else ''}{sample_name}_{tRNA_annotation}.bam"
                bam_filepath = os.path.join(bam_dir, bam_filename)

                # Check if BAM file exists and overwrite is False
                if os.path.exists(bam_filepath) and not overwrite:
                    print(f"BAM file already exists: {bam_filepath}. Skipping...")
                    continue

                # Open the input JSON file
                try:
                    with bz2.open(json_file, 'rt', encoding="utf-8") as json_fh:
                        alignments = json_stream.load(json_fh)
                        alignments = alignments.persistent()

                        # Filter reads based on tRNA_annotation
                        filtered_reads = {
                            read_id: align_data
                            for read_id, align_data in alignments.items()
                            if align_data['name'] == tRNA_annotation and align_data['aligned']
                        }

                    # Write to BAM
                    with AlignmentFile(bam_filepath, "wb", header={"HD": {"VN": "1.0"}, "SQ": [{"LN": 1, "SN": tRNA_annotation}]}) as bam_out:
                        for read_id, align_data in filtered_reads.items():
                            bam_out.write({
                                'qname': read_id,
                                'flag': 0,
                                'rname': tRNA_annotation,
                                'pos': align_data['dpos'][0],
                                'mapq': align_data['score'],
                                'cigar': [],  # Populate with appropriate CIGAR string
                                'rnext': '*',
                                'pnext': 0,
                                'tlen': 0,
                                'seq': align_data['qseq'],
                                'qual': '*'  # Quality score
                            })

                    print(f"BAM file created: {bam_filepath}")

                except Exception as e:
                    print(f"Error processing {sample_name}: {e}")

def row_exist_or_none(row, col):
    try:
        return(row[col])
    except:
        return(None)

def UMI_exp(c, cUMI, N_bins):
    b = 4**9*2
    N_umi_exp = N_bins*(1-((N_bins-1) / N_bins)**c)
    return(100 * cUMI / N_umi_exp)

