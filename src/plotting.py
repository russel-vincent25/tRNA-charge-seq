import os, shutil, bz2, warnings, copy, contextlib
import Bio.Data.CodonTable
import pandas as pd
import numpy as np
from mpire import WorkerPool
from Bio import Align

### Plotting imports and settings ###
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.patches import StepPatch
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import logomaker as lm
palette = list(mcolors.TABLEAU_COLORS.keys())
sns.set_theme(style="ticks", palette="muted")
sns.set_context("talk")



class TRNA_plot:
    '''
    Class to generate some standard plots
    from the tRNAseq data such as charge per codon,
    read per million (RPM) etc.

    Keyword arguments:
    sample_df -- sample_df pandas dataframe. Can also be pulled from default location (default None)
    pull_default -- Pull sample_df from default location in the alignment folder
    stats_fnam -- Filename for the aggregated stats csv file made during stats collection. If None is specified trying to look it up in the stats collection dir under ALL_stats_aggregate.csv (default None)
    overwrite_dir -- Overwrite old plotting folder if any exists (default False)
    charge_count -- Column to use for charge calculation. Either count or UMIcount (default 'count')
    RPM_count -- Column to use for RPM and coverage calculations. Either count or UMIcount (default 'UMIcount')
    excl_align_gap -- Exclude all alignments with a gap (default False)
    excl_09_fmax -- Exclude all alignments with an alignment score less than 90% of maximum for the given length alignment (default False)
    '''
    def __init__(self, dir_dict, sample_df=None, stats_fnam=None, \
                 pull_default=False, overwrite_dir=False, \
                 charge_count='count', RPM_count='UMIcount', \
                 excl_align_gap=False, excl_09_fmax=False):
        self.stats_csv_header = ['readID', 'common_seq', 'sample_name_unique', \
                                 'sample_name', 'replicate', 'barcode', 'species', 'tRNA_annotation', \
                                 'align_score', 'fmax_score', 'Ndeletions', 'Ninsertions', \
                                 'unique_annotation', 'tRNA_annotation_len', \
                                 'align_5p_idx', 'align_3p_idx', 'align_5p_nt', 'align_3p_nts', \
                                 'codon', 'anticodon', 'amino_acid', '5p_cover', '3p_cover', \
                                 '5p_non-temp', '3p_non-temp', '5p_UMI', '3p_BC', \
                                 'align_gap', 'fmax_score>0.9', 'UMIcount', 'count']
        # Adding new classification fields
        self.stats_csv_header.extend([
            "is_5p_tRF", "is_degraded_tRNA", "is_pre_tRNA",
            "5p_non_temp_seq", "3p_non_temp_seq"
        ])
        self.stats_csv_header_type = [str, bool, str, str, int, str, str, str, int, float, \
                                      int, int, str, int, int, int, str, str, str, str, str, \
                                      bool, bool, str, str, str, str, bool, bool, int, int]
        # Adding new classification fields
        self.stats_csv_header_type.extend([
            bool, bool, bool, str, str
        ])
        self.stats_csv_header_td = {nam:tp for nam, tp in zip(self.stats_csv_header, self.stats_csv_header_type)}
        
        self.stats_agg_cols = ['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'species', \
                               'tRNA_annotation', 'tRNA_annotation_len', 'unique_annotation', \
                               '5p_cover','3p_cover','align_5p_idx', 'align_3p_idx', 'align_3p_nts','3p_non-temp','codon', 'anticodon', 'amino_acid', \
                               'align_gap', 'fmax_score>0.9', 'UMIcount', 'count', 'UMI_percent_exp']
        # Adding the 3' end charge categories and new classifications
        self.stats_agg_cols.extend([
            "Total_CA", "Total_GA", "Total_CC", "Total_CG",
            "Total_5p_tRFs", "Total_Degraded", "Total_Pre_tRNA"
        ])            
        self.stats_agg_cols_type = [str, str, int, str, str, str, int, bool, \
                                    bool, bool, int, int, str, str , str, str, str, bool, bool, int, int, float]
        self.stats_agg_cols_type.extend([
            int, int, int, int, int, int, int
        ])
        self.stats_agg_cols_td = {nam:tp for nam, tp in zip(self.stats_agg_cols, self.stats_agg_cols_type)}

        # Input:
        self.dir_dict = dir_dict
        self.charge_filt = dict()
        self.excl_align_gap = excl_align_gap
        self.excl_09_fmax = excl_09_fmax
        if charge_count in ['count', 'UMIcount'] and RPM_count in ['count', 'UMIcount']:
            self.charge_count_col = charge_count
            self.RPM_count_col = RPM_count
        else:
            raise Exception('"charge_count" and "RPM_count" input must be either "count" or "UMIcount"')

        self.stats_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['stats_dir'])
        # Attempt to load sample_df from default path:
        if pull_default:
            align_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['align_dir'])
            sample_df_path = '{}/{}'.format(align_dir_abs, 'sample_stats.xlsx')
            if os.path.isfile(sample_df_path):
                sample_df = pd.read_excel(sample_df_path, index_col=0)
            else:
                raise Exception('Default "sample_df" could not be found: {}'.format(sample_df_path))
        elif sample_df is None:
            raise Exception('A "sample_df" must be specified or "pulled from default".')
        self.sample_df = sample_df

        # Load aggregated CSV data:
        if stats_fnam is None:
            stats_fnam = '{}/ALL_stats_aggregate.csv'.format(self.stats_dir_abs)
        self.all_stats = pd.read_csv(stats_fnam, keep_default_na=False, dtype=self.stats_agg_cols_td)
        # Get rid of 1/2 in mito Leu amino acid name:
        self.all_stats['amino_acid'] = [AA[:-1] if AA[-1]=='1' or AA[-1]=='2' else AA for AA in self.all_stats['amino_acid'].values]
        # Translate codon into single letter amino acid:
        self.all_stats['AA_letter'] = [AAA2A[AAA] for AAA in self.all_stats['amino_acid'].values]
        # Add amino acid - codon string:
        self.all_stats['AA_codon'] = [AA + '-' + codon for codon, AA in zip(self.all_stats['codon'].values, self.all_stats['amino_acid'].values)]
        # tRNA annotation short form:
        tshort_lst = list()
        for an in self.all_stats['tRNA_annotation'].values:
            an_short = list()
            for an_p in an.split('@'):
                if '_tRX' in an_p:
                    an_short.append('-'.join(an_p.split('-')[1:]) + '-X')
                else:
                    an_short.append('-'.join(an_p.split('-')[1:]))
            tshort_lst.append('@'.join(an_short))
        self.all_stats['tRNA_anno_short'] = tshort_lst

        # Create single codon filter, to filter out sequences that map to
        # tRNA sequences with different codon/anticodon:
        single_codon = list()
        for anno_str, anticodon in zip(self.all_stats['tRNA_annotation'].values, self.all_stats['anticodon'].values):
            sc = True
            anno_list = anno_str.split('@')
            for anno in anno_list:
                if not anno.split('-')[2] == anticodon:
                    sc = False
            single_codon.append(sc)
        self.all_stats['single_codon'] = single_codon
        # Create single amino acid filter, to filter out sequences that map to
        # tRNA sequences with different codon/anticodon:
        single_aa = list()
        for anno_str, amino_acid in zip(self.all_stats['tRNA_annotation'].values, self.all_stats['amino_acid'].values):
            sa = True
            anno_list = anno_str.split('@')
            for anno in anno_list:
                # Use "in" instead of "==" to accommodate Leu1/2 and Ser1/2:
                if not amino_acid in anno.split('-')[1]:
                    sa = False
            single_aa.append(sa)
        self.all_stats['single_aa'] = single_aa
        self.all_stats['mito_codon'] = ['mito_tRNA' in anno for anno in self.all_stats['tRNA_annotation'].values]
        self.all_stats['Syn_ctr'] = ['Synthetic' in anno and sp != 'ecoli' for sp, anno in zip(self.all_stats['species'].values, self.all_stats['tRNA_annotation'].values)]
        # Add new columns to stats_agg:
        self.stats_agg_cols = ['sample_name_unique', 'sample_name', 'replicate', 'barcode', 'species', \
                               'tRNA_annotation', 'tRNA_anno_short', 'tRNA_annotation_len', 'unique_annotation', \
                               '5p_cover', '3p_cover','align_3p_nts', 'codon', 'anticodon', 'amino_acid', 'AA_letter', \
                               'AA_codon', 'single_codon', 'single_aa', 'mito_codon', 'Syn_ctr', \
                               'align_gap', 'fmax_score>0.9', 'UMIcount', 'count','is_5p_tRF', \
                               'is_degraded_tRNA', 'is_pre_tRNA','5p_non_temp_seq', '3p_non_temp_seq']
        self.stats_agg_cols_td = {'sample_name_unique': str, 'sample_name': str, 'replicate': int, \
                                  'barcode': str, 'species': str, 'tRNA_annotation': str, 'tRNA_anno_short': str, \
                                  'tRNA_annotation_len': int, 'unique_annotation': bool, \
                                  '5p_cover': bool, '3p_cover':bool, 'align_3p_nts': str, 'codon': str, 'anticodon': str, 'amino_acid': str, 'AA_letter': str, \
                                  'AA_codon': str, 'single_codon': bool, 'single_aa': bool, 'mito_codon': bool, 'Syn_ctr': bool, \
                                  'align_gap': bool, 'fmax_score>0.9': bool, 'UMIcount': int, 'count': int, \
                                  'is_5p_tRF': bool, 'is_degraded_tRNA': bool, 'is_pre_tRNA': bool, \
                                  '5p_non_temp_seq': str, '3p_non_temp_seq': str}
        
        # Reorder columns:
        self.all_stats = self.all_stats.loc[:, self.stats_agg_cols].copy()
        # Calculate charge and create dataframe:
        self._get_charge_df()

        # Make output folder:
        self._make_dir(overwrite_dir)

    def _make_dir(self, overwrite):
        # Create folder for files:
        self.plotting_dir_abs = '{}/{}/{}'.format(self.dir_dict['NBdir'], self.dir_dict['data_dir'], self.dir_dict['plotting_dir'])
        try:
            os.mkdir(self.plotting_dir_abs)
        except:
            if overwrite:
                shutil.rmtree(self.plotting_dir_abs)
                os.mkdir(self.plotting_dir_abs)
            else:
                print('Folder exists and overwrite set to false... Doing nothing.')

    def _get_charge_df(self):
        # Filter and group to only extract rows
        # valid for charge determination:
        row_mask = (self.all_stats['count'] > 0) # all True (dummy)
        if self.excl_align_gap:
            row_mask &= (~self.all_stats['align_gap'])
        if self.excl_09_fmax:
            row_mask &= (self.all_stats['fmax_score>0.9'])
        charge_df = self.all_stats.loc[row_mask, self.stats_agg_cols].copy()
        # Rearrange the columns of charge_df such that count and UMIcount are last:
        charge_df = charge_df[[col for col in charge_df.columns if col not in ['count', 'UMIcount']] + ['count', 'UMIcount']]
        # Group by all columns expcept count and UMIcount:
        charge_df = charge_df.groupby(self.stats_agg_cols[:-2], as_index=False).agg({'count': "sum", 'UMIcount': "sum"}).reset_index(drop=True)

        # Count A, C, G endings:
        charge_df['CA_count'] = [ct if nt == 'CA' else 0 for ct, nt in zip(charge_df[self.charge_count_col], charge_df['align_3p_nts'])]
        charge_df['CC_count'] = [ct if nt == 'CC' else 0 for ct, nt in zip(charge_df[self.charge_count_col], charge_df['align_3p_nts'])]
        charge_df['GA_count'] = [ct if nt == 'GA' else 0 for ct, nt in zip(charge_df[self.charge_count_col], charge_df['align_3p_nts'])]
        charge_df['CG_count'] = [ct if nt == 'CG' else 0 for ct, nt in zip(charge_df[self.charge_count_col], charge_df['align_3p_nts'])]
        # Count 5' tRFs, degraded and pre-tRNAs for boolean columns:
        charge_df['5p_tRF'] = [ct if tRF else 0 for ct, tRF in zip(charge_df[self.charge_count_col], charge_df['is_5p_tRF'])]
        charge_df['degraded'] = [ct if deg else 0 for ct, deg in zip(charge_df[self.charge_count_col], charge_df['is_degraded_tRNA'])]
        charge_df['pre_tRNA'] = [ct if pre else 0 for ct, pre in zip(charge_df[self.charge_count_col], charge_df['is_pre_tRNA'])]
        
        # Group transcripts with different 3p nt.
        # then calculate charge:
        charge_df_cols = copy.deepcopy(self.stats_agg_cols)
        charge_df_cols.remove('align_3p_nts')
        charge_df_cols.remove('align_gap')
        charge_df_cols.remove('fmax_score>0.9')
        charge_df_cols.extend(['CA_count', 'CC_count', 'GA_count', 'CG_count', '5p_tRF', 'degraded', 'pre_tRNA'])

        charge_df = charge_df.groupby(charge_df_cols[:-9], as_index=False).agg({'count': "sum", \
                                                                                'UMIcount': "sum", \
                                                                                'CA_count': "sum", \
                                                                                'CC_count': "sum", \
                                                                                'GA_count': "sum", \
                                                                                'CG_count': "sum", \
                                                                                '5p_tRF': "sum", \
                                                                                'degraded': "sum", \
                                                                                'pre_tRNA': "sum"}).reset_index(drop=True)
        # Safe division with NaN for zero denominators
        charge_df['charge_canonical'] = charge_df.apply(
            lambda row: 100 * row['CA_count'] / (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded'])
            if (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded']) > 0 else np.nan, axis=1)
        charge_df['charge_non-canonical'] = charge_df.apply(
            lambda row: 100 * row['GA_count'] / (row['GA_count'] + row['CG_count'] + row['5p_tRF'] + row['degraded'])
            if (row['GA_count'] + row['CG_count'] + row['5p_tRF'] + row['degraded']) > 0 else np.nan, axis=1)
        charge_df['tRNA_fragment'] = charge_df.apply(
            lambda row: (row['5p_tRF'] + row['degraded']) / row['count'] if row['count'] > 0 else np.nan, axis=1)

        # Add the sample total count to the rows:
        df_count = charge_df[~charge_df['Syn_ctr']].groupby(['sample_name_unique'], as_index=False).agg({self.RPM_count_col: "sum"}).reset_index(drop=True)
        charge_df = charge_df.merge(df_count, on='sample_name_unique', suffixes=('', '_sample_tot'))
        # Calculated the RPM and get rid of the total count:
        charge_df['RPM'] = charge_df[self.RPM_count_col] / (charge_df['{}_sample_tot'.format(self.RPM_count_col)] / 1e6)
        charge_df = charge_df.drop(columns=['{}_sample_tot'.format(self.RPM_count_col)])

        # Filter and group by same amino acid:
        aa_mask = charge_df['single_aa']
        charge_df_aa = charge_df[aa_mask].groupby(['sample_name_unique', 'sample_name', 'replicate', \
                                                   'barcode', 'amino_acid', 'AA_letter', 'mito_codon', \
                                                   'Syn_ctr'], as_index=False).agg({'count': "sum", \
                                                                                    'UMIcount': "sum", \
                                                                                    'CA_count': "sum", \
                                                                                    'CC_count': "sum", \
                                                                                    'GA_count': "sum", \
                                                                                    'CG_count': "sum", \
                                                                                    '5p_tRF': "sum", \
                                                                                    'degraded': "sum", \
                                                                                    'pre_tRNA': "sum", \
                                                                                    'RPM': "sum"}).reset_index(drop=True)
        # Safe division with NaN for zero denominators
        charge_df_aa['charge_canonical'] = charge_df_aa.apply(
            lambda row: 100 * row['CA_count'] / (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded'])
            if (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded']) > 0 else np.nan, axis=1)
        charge_df_aa['charge_non-canonical'] = charge_df_aa.apply(
            lambda row: 100 * row['GA_count'] / (row['GA_count'] + row['CG_count'] + row['5p_tRF'] + row['degraded'])
            if (row['GA_count'] + row['CG_count'] + row['5p_tRF'] + row['degraded']) > 0 else np.nan, axis=1)
        charge_df_aa['tRNA_fragment'] = charge_df_aa.apply(
            lambda row: (row['5p_tRF'] + row['degraded']) / row['count'] if row['count'] > 0 else np.nan, axis=1)
        self.charge_filt['aa'] = charge_df_aa

        # Filter and group by same codon:
        cd_mask = charge_df['single_codon']
        charge_df_cd = charge_df[cd_mask].groupby(['sample_name_unique', 'sample_name', 'replicate', \
                                                   'barcode', 'codon', 'anticodon', 'AA_codon', \
                                                   'amino_acid', 'AA_letter', 'mito_codon', \
                                                   'Syn_ctr'], as_index=False).agg({'count': "sum", \
                                                                                    'UMIcount': "sum", \
                                                                                    'CA_count': "sum", \
                                                                                    'CC_count': "sum", \
                                                                                    'GA_count': "sum", \
                                                                                    'CG_count': "sum", \
                                                                                    '5p_tRF': "sum", \
                                                                                    'degraded': "sum", \
                                                                                    'pre_tRNA': "sum", \
                                                                                    'RPM': "sum"}).reset_index(drop=True)
        # Safe division with NaN for zero denominators
        charge_df_cd['charge_canonical'] = charge_df_cd.apply(
            lambda row: 100 * row['CA_count'] / (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded'])
            if (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded']) > 0 else np.nan, axis=1)
        charge_df_cd['charge_non-canonical'] = charge_df_cd.apply(
            lambda row: 100 * row['GA_count'] / (row['GA_count'] + row['CG_count'] + row['5p_tRF'] + row['degraded'])
            if (row['GA_count'] + row['CG_count'] + row['5p_tRF'] + row['degraded']) > 0 else np.nan, axis=1)
        charge_df_cd['tRNA_fragment'] = charge_df_cd.apply(
            lambda row: (row['5p_tRF'] + row['degraded']) / row['count'] if row['count'] > 0 else np.nan, axis=1)
        self.charge_filt['codon'] = charge_df_cd

        # Filter and group by same transcript:
        tr_mask = charge_df['unique_annotation']
        charge_df_tr = charge_df[tr_mask].groupby(['sample_name_unique', 'sample_name', 'replicate', \
                                                   'barcode', 'tRNA_annotation', 'tRNA_anno_short', \
                                                   'tRNA_annotation_len', 'codon', 'anticodon', 'AA_codon', \
                                                   'amino_acid', 'AA_letter', 'mito_codon', \
                                                   'Syn_ctr'], as_index=False).agg({'count': "sum", \
                                                                                    'UMIcount': "sum", \
                                                                                    'CA_count': "sum", \
                                                                                    'CC_count': "sum", \
                                                                                    'GA_count': "sum", \
                                                                                    'CG_count': "sum", \
                                                                                    '5p_tRF': "sum", \
                                                                                    'degraded': "sum", \
                                                                                    'pre_tRNA': "sum", \
                                                                                    'RPM': "sum"}).reset_index(drop=True)

        # Safe division with NaN for zero denominators
        charge_df_tr['charge_canonical'] = charge_df_tr.apply(
            lambda row: 100 * row['CA_count'] / (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded'])
            if (row['CA_count'] + row['CC_count'] + row['5p_tRF'] + row['degraded']) > 0 else np.nan, axis=1)
        charge_df_tr['charge_non-canonical'] = charge_df_tr.apply(
            lambda row: 100 * row['GA_count'] / (row['GA_count'] + row['CG_count'] + row['5p_tRF'] + row['degraded'])
            if (row['GA_count'] + row['CG_count'] + row['5p_tRF'] + row['degraded']) > 0 else np.nan, axis=1)
        charge_df_tr['tRNA_fragment'] = charge_df_tr.apply(
            lambda row: (row['5p_tRF'] + row['degraded']) / row['count'] if row['count'] > 0 else np.nan, axis=1)
        self.charge_filt['tr'] = charge_df_tr
        self.charge_df = charge_df

    def write_charge_df(self, df_type='transcript', fnam='charge_df'):
        '''qwerty'''
        if df_type == 'aa':
            charge_df_type = self.charge_filt['aa'].copy()
        elif df_type == 'codon':
            charge_df_type = self.charge_filt['codon'].copy()
        elif df_type == 'transcript':
            charge_df_type = self.charge_filt['tr'].copy()
        else:
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "aa", "codon" or "transcript".'.format(plot_type))
        fnam_abs = '{}/{}.csv'.format(self.plotting_dir_abs, fnam)
        charge_df_type.to_csv(fnam_abs, header=True, index=False)

    def plot_Syn_ctr(self, plot_name='ecoli-ctr_charge_plot', sample_list=None,
                       charge_plot=True, min_obs=1, \
                       sample_list_exl=None, bc_list_exl=None):

        if charge_plot:
            y_axis = 'charge'
            count_col = self.charge_count_col
        else:
            y_axis = 'RPM'
            count_col = self.RPM_count_col

        if sample_list is None:
            sample_list = list(set(self.sample_df['sample_name_unique']))

        # Sample rows selected:
        charge_df_type = self.charge_filt['tr'].copy()
        sample_mask = np.array([False]*len(charge_df_type))
        for unam in sample_list:
            sample_mask |= (charge_df_type['sample_name_unique'] == unam)

        # Exclude unique samples or barcodes:
        if not sample_list_exl is None:
            for unam in sample_list_exl:
                sample_mask &= (charge_df_type['sample_name_unique'] != unam)
        if not bc_list_exl is None:
            for bc in bc_list_exl:
                sample_mask &= (charge_df_type['barcode'] != bc)

        # Only take rows with Ecoli controls:
        sample_mask &= charge_df_type['Syn_ctr'] & (charge_df_type[count_col] >= min_obs)
        charge_sample = charge_df_type[sample_mask].copy()
        
        # Different settings for different plot types:
        Ngrp = len(set(charge_sample['sample_name_unique']))
        x_list = sorted(set(charge_sample['sample_name_unique'].values), key=str.casefold)

        # Make the plot:
        fig, ax = plt.subplots(1, 1, figsize=(0.6*Ngrp, 6))
        g1 = sns.barplot(ax=ax, x='sample_name_unique', y=y_axis, order=x_list, hue='tRNA_anno_short', \
                         data=charge_sample, capsize=.1, errwidth=2, \
                         edgecolor='black', linewidth=2, alpha=0.85)

        # Set axis/title text:
        g1.set_title('E.coli control samples')
        g1.set_xticklabels(g1.get_xticklabels(), rotation=90)
        g1.grid(True, axis='y')
        g1.set_xlabel('');
        
        if charge_plot is True:
            g1.set_ylabel('Charge (%)');
            g1.set(ylim=[0, 100])
        else:
            g1.set_ylabel('Abundance (RPM)');

        # Write to PDF file:
        plot_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        fig.tight_layout()
        fig.savefig(plot_fnam, bbox_inches='tight')
        plt.close(fig)

    def plot_abundance(self, plot_type='aa', group=False, charge_plot=False, \
        plot_name='abundance_plot_aa', min_obs=100, sample_list=None, verbose=True, \
        sample_list_exl=None, bc_list_exl=None, tr_short=True):

        if plot_type == 'aa':
            charge_df_type = self.charge_filt['aa'].copy()
        elif plot_type == 'codon':
            charge_df_type = self.charge_filt['codon'].copy()
        elif plot_type == 'transcript':
            charge_df_type = self.charge_filt['tr'].copy()
        else:
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "aa", "codon" or "transcript".'.format(plot_type))
        
        if charge_plot:
            y_axis = 'charge'
            count_col = self.charge_count_col
        else:
            y_axis = 'RPM'
            count_col = self.RPM_count_col

        if sample_list is None:
            sample_list = list(self.sample_df['sample_name'].drop_duplicates())

        if verbose:
            print('\nNow plotting sample/group:', end='')

        # Use a 40 color colormap:
        cmap_b = mpl.colormaps['tab20']
        cmap_c = mpl.colormaps['tab20']
        cmap_b.colors = cmap_b.colors + cmap_c.colors
        cmap_b.N = 40
        cmap_b._i_bad = 42
        cmap_b._i_over = 40
        cmap_b._i_under = 0

        snam_df = self.sample_df.loc[:, ['sample_name', 'plot_group', 'hue_name', 'hue_value', 'hue_order']].drop_duplicates()
        charge_df_type = charge_df_type.merge(snam_df, on='sample_name')
        if group:
            # Find all sample groups:
            plot_group_set = [set(pg.split(' @&@ ')) for pg in self.sample_df.loc[:, 'plot_group'].values]
            plot_iter = sorted({list(pg)[0] for pg in plot_group_set if len(pg) == 1})
        else:
            plot_iter = snam_df['sample_name'].drop_duplicates()

        # Print each plot to the same PDF file:
        charge_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(charge_fnam) as pp:
            # Loop through and generate plots for each sample:
            for sg_name in plot_iter:
                if group is False and not sg_name in sample_list:
                    continue
                print('  {}'.format(sg_name), end='')

                # Sample rows selected:
                if group:
                    # Pick out samples names in the group:
                    snam_mask = np.array([sg_name in pg for pg in plot_group_set])
                    snam_set = set(self.sample_df.loc[snam_mask, 'sample_name'])
                    sample_mask = (charge_df_type['sample_name'].isin(snam_set)) & (charge_df_type[count_col] >= min_obs)
                else:
                    sample_mask = (charge_df_type['sample_name'] == sg_name) & (charge_df_type[count_col] >= min_obs)

                # Exclude unique samples or barcodes:
                if not sample_list_exl is None:
                    for unam in sample_list_exl:
                        sample_mask &= (charge_df_type['sample_name_unique'] != unam)
                if not bc_list_exl is None:
                    for bc in bc_list_exl:
                        sample_mask &= (charge_df_type['barcode'] != bc)

                charge_sample = charge_df_type[sample_mask].copy()
                # Plot separate for mito/cyto:
                mask_cyto = (~charge_sample['Syn_ctr']) & (~charge_sample['mito_codon'])
                mask_mito = (~charge_sample['Syn_ctr']) & (charge_sample['mito_codon'])

                # Different settings for different plot types:
                Ngrp = len(set(charge_sample['sample_name']))
                hue_order = [t[0] for t in sorted(dict(zip(charge_sample['hue_value'].values, charge_sample['hue_order'].values)).items(), key=lambda x: x[1])]
                if plot_type == 'aa':
                    x_axis = 'amino_acid'
                    x_list_cyto = sorted(set(charge_sample['amino_acid'][mask_cyto].values), key=str.casefold)
                    x_list_mito = sorted(set(charge_sample['amino_acid'][mask_mito].values), key=str.casefold)
                    AA_set = set(x_list_cyto+x_list_mito)
                    AAi = {aa:i for i, aa in enumerate(sorted(AA_set))}
                    colors_cyto = {c: cmap_b(AAi[c]) for c in x_list_cyto}
                    colors_mito = {c: cmap_b(AAi[c]) for c in x_list_mito}
                    if sum(mask_mito) > 0:
                        fig = plt.figure(figsize=(8*Ngrp, 9))
                        gs = fig.add_gridspec(2, 4)
                        ax1 = fig.add_subplot(gs[0, :])
                        ax2 = fig.add_subplot(gs[1, :])
                    else:
                        fig = plt.figure(figsize=(8*Ngrp, 4.5))
                        gs = fig.add_gridspec(1, 4)
                        ax1 = fig.add_subplot(gs[0, :])
                elif plot_type == 'codon':
                    x_axis = 'AA_codon'
                    x_list_cyto = sorted(set(charge_sample['AA_codon'][mask_cyto].values), key=str.casefold)
                    x_list_mito = sorted(set(charge_sample['AA_codon'][mask_mito].values), key=str.casefold)
                    AA_set = {c.split('-')[0] for c in x_list_cyto+x_list_mito}
                    AAi = {aa:i for i, aa in enumerate(sorted(AA_set))}
                    colors_cyto = {c: cmap_b(AAi[c.split('-')[0]]) for c in x_list_cyto}
                    colors_mito = {c: cmap_b(AAi[c.split('-')[0]]) for c in x_list_mito}
                    if sum(mask_mito) > 0:
                        fig = plt.figure(figsize=(16*Ngrp, 9))
                        gs = fig.add_gridspec(2, 4)
                        ax1 = fig.add_subplot(gs[0, :])
                        ax2 = fig.add_subplot(gs[1, 1:3])
                    else:
                        fig = plt.figure(figsize=(16*Ngrp, 4.5))
                        gs = fig.add_gridspec(1, 4)
                        ax1 = fig.add_subplot(gs[0, :])
                elif plot_type == 'transcript':
                    if tr_short:
                        x_axis = 'tRNA_anno_short'
                    else:
                        x_axis = 'tRNA_annotation'
                    mask_cyto = (~charge_sample['Syn_ctr']) & (~charge_sample['mito_codon'])
                    mask_mito = (~charge_sample['Syn_ctr']) & (charge_sample['mito_codon'])
                    x_list_cyto = sorted(set(charge_sample[x_axis][mask_cyto].values), key=str.casefold)
                    x_list_mito = sorted(set(charge_sample[x_axis][mask_mito].values), key=str.casefold)
                    colors_cyto = None
                    colors_mito = None
                    if sum(mask_mito) > 0:
                        fig = plt.figure(figsize=(45*Ngrp, 15))
                        gs = fig.add_gridspec(2, 10)
                        ax1 = fig.add_subplot(gs[0, :])
                        ax2 = fig.add_subplot(gs[1, 4:6])
                    else:
                        fig = plt.figure(figsize=(45*Ngrp, 7.5))
                        gs = fig.add_gridspec(1, 10)
                        ax1 = fig.add_subplot(gs[0, :])

                # Plot either as grouped or single bars:
                if group:
                    # Cyto tRNAs
                    g1 = sns.barplot(ax=ax1, x=x_axis, y=y_axis, hue='hue_value', hue_order=hue_order, order=x_list_cyto, data=charge_sample[mask_cyto], capsize=.025, errwidth=2, edgecolor='black', linewidth=2, alpha=0.8)
                    if sum(mask_mito) > 0:
                        # Mito tRNAs
                        g2 = sns.barplot(ax=ax2, x=x_axis, y=y_axis, hue='hue_value', hue_order=hue_order, order=x_list_mito, data=charge_sample[mask_mito], capsize=.025, errwidth=2, edgecolor='black', linewidth=2, alpha=0.8)
                        # Legend:
                        g2.legend_.remove()
                    old_legend = g1.legend_
                    handles = old_legend.legendHandles
                    labels = hue_order
                    title = charge_sample['hue_name'].values[0]
                    g1.legend(handles, labels, title=title, bbox_to_anchor=(1.01,1), borderaxespad=0)
                else:
                    # Cyto tRNAs
                    g1 = sns.barplot(ax=ax1, x=x_axis, y=y_axis, order=x_list_cyto, data=charge_sample[mask_cyto], capsize=.1, errwidth=2, edgecolor='black', linewidth=2, alpha=0.85, palette=colors_cyto)
                    if sum(mask_mito) > 0:
                        # Mito tRNAs
                        g2 = sns.barplot(ax=ax2, x=x_axis, y=y_axis, order=x_list_mito, data=charge_sample[mask_mito], capsize=.1, errwidth=2, edgecolor='black', linewidth=2, alpha=0.85, palette=colors_mito)

                # Set axis/title text:
                g1.set_title('Cytoplasmic tRNA for sample/group {}'.format(sg_name))
                g1.set_xticklabels(g1.get_xticklabels(), rotation=90)
                g1.grid(True, axis='y')
                g1.set_xlabel('');
                if sum(mask_mito) > 0:
                    g2.set_title('Mitochondrial tRNA for sample/group {}'.format(sg_name))
                    g2.set_xticklabels(g2.get_xticklabels(), rotation=90)
                    g2.grid(True, axis='y')
                    g2.set_xlabel('');

                if charge_plot is True:
                    g1.set_ylabel('Charge (%)');
                    g1.set(ylim=[0, 100])
                    if sum(mask_mito) > 0:
                        g2.set_ylabel('Charge (%)');
                        g2.set(ylim=[0, 100])
                else:
                    g1.set_ylabel('Abundance (RPM)');
                    if sum(mask_mito) > 0:
                        g2.set_ylabel('Abundance (RPM)');

                # Write to PDF file:
                fig.tight_layout()
                pp.savefig(fig, bbox_inches='tight')
                plt.close(fig)

    def plot_abundance_corr(self, \
                            sample_pairs=None, \
                            sample_pairs_col='sample_name', \
                            sample_unique_pairs=None, \
                            plot_type='aa', charge_plot=False, \
                            log=False, plot_name='abundance_corr_plot_aa', \
                            min_obs=100, sample_list=None, verbose=True, \
                            sample_list_exl=None, \
                            bc_list_exl=None, one2one_corr=False):

        if plot_type == 'aa':
            charge_df_type = self.charge_filt['aa'].copy()
            merge_on = ['amino_acid', 'mito_codon', 'Syn_ctr']
        elif plot_type == 'codon':
            charge_df_type = self.charge_filt['codon'].copy()
            merge_on = ['AA_codon', 'mito_codon', 'Syn_ctr']
        elif plot_type == 'transcript':
            charge_df_type = self.charge_filt['tr'].copy()
            merge_on = ['tRNA_annotation', 'mito_codon', 'Syn_ctr']
        else:
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "aa", "codon" or "transcript".'.format(plot_type))

        if charge_plot:
            metric = 'charge'
            count_col = self.charge_count_col
        else:
            metric = 'RPM'
            count_col = self.RPM_count_col

        # Enforce minimum observations,
        # and not Ecoli control:
        min_obs_mask = (charge_df_type[count_col] > min_obs) & ~charge_df_type['Syn_ctr']
        charge_df_type = charge_df_type[min_obs_mask]

        # Handle if input is unique sample names
        # or sample names with replicates to merge:
        if not sample_pairs is None:
            assert(len(sample_pairs[0]) == len(sample_pairs[1]))
            # Exclude samples/barcodes:
            mask_exl = np.array([True]*len(self.sample_df))
            if not sample_list_exl is None:
                for snam in sample_list_exl:
                    mask_exl &= (self.sample_df['sample_name_unique'] != snam)
            if not bc_list_exl is None:
                for bc in bc_list_exl:
                    mask_exl &= (self.sample_df['barcode'] != bc)

            # Convert sample names to lists of unique sample names:
            pair_list = [[], []]
            name_list = [[], []]
            for sp1, sp2 in zip(*sample_pairs):
                mask1 = mask_exl & (self.sample_df[sample_pairs_col] == sp1)
                mask2 = mask_exl & (self.sample_df[sample_pairs_col] == sp2)
                if mask1.sum() > 0 and mask2.sum() > 0:
                    pair_list[0].append(list(self.sample_df.loc[mask1, 'sample_name_unique'].values))
                    pair_list[1].append(list(self.sample_df.loc[mask2, 'sample_name_unique'].values))
                    name_list[0].append(sp1)
                    name_list[1].append(sp2)

        elif not sample_unique_pairs is None:
            assert(len(sample_unique_pairs[0]) == len(sample_unique_pairs[1]))
            pair_list = [[], []]
            name_list = [[], []]
            for sp1, sp2 in zip(*sample_unique_pairs):
                pair_list[0].append([sp1])
                pair_list[1].append([sp2])
                name_list[0].append(sp1)
                name_list[1].append(sp2)
        else:
            raise Exception('Neither sample_pairs nor sample_unique_pairs lists were provided.')

        if sample_list is None:
            sample_list = list(self.sample_df['sample_name'].drop_duplicates())

        if verbose:
            print('\nNow plotting sample pairs:', end='')

        # Print each plot to the same PDF file:
        corr_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(corr_fnam) as pp:
            # Plot each pair:
            for name_idx, sp1l_sp2l in enumerate(zip(*pair_list)):
                sp1l, sp2l = sp1l_sp2l
                nam1, nam2 = name_list[0][name_idx], name_list[1][name_idx]
                # Select rows for all samples in each group:
                mask1 = charge_df_type['sample_name_unique'].isin(sp1l)
                mask2 = charge_df_type['sample_name_unique'].isin(sp2l)

                if mask1.sum() == 0 or mask2.sum() == 0:
                    continue
                print('  ({} - {})'.format(nam1, nam2), end='')

                # If charge is plotted both samples
                # have to be present, if RPM is plotted
                # it can be zero in one sample and >0 in another:
                if not charge_plot and min_obs == 0:
                    merge_how = 'outer'
                else:
                    merge_how = 'inner'

                # Merge all samples in each group:
                sele_cols = merge_on + ['RPM', 'charge']
                m1 = charge_df_type[mask1].copy()
                m1 = m1[sele_cols].groupby(sele_cols[:-2], as_index=False).agg({'RPM': 'mean', 'charge': 'mean'}).reset_index(drop=True)
                m2 = charge_df_type[mask2].copy()
                m2 = m2[sele_cols].groupby(sele_cols[:-2], as_index=False).agg({'RPM': 'mean', 'charge': 'mean'}).reset_index(drop=True)
                
                # Merge the two:
                m12 = pd.merge(
                    m1,
                    m2,
                    how=merge_how,
                    on=merge_on,
                    left_index=False,
                    right_index=False,
                    sort=True,
                    suffixes=("_1", "_2"),
                    copy=True,
                    indicator=False,
                    validate=None,
                )
                m12 = m12.fillna(0)

                # Handle if log transform is requested:
                if log:
                    mask0 = (m12[metric+'_1'] > 0) & (m12[metric+'_2'] > 0)
                    m12 = m12[mask0]
                    m12[metric+'_1'] = np.log10(m12[metric+'_1'])
                    m12[metric+'_2'] = np.log10(m12[metric+'_2'])

                # Plot the scatterplot:
                fig, ax = plt.subplots(1, 1, figsize=(7, 6))
                g1 = sns.regplot(ax=ax, data=m12, x=metric+'_1', y=metric+'_2')
                if not charge_plot or one2one_corr:
                    xy_min = [min(g1.get_xlim()), min(g1.get_ylim())]
                    xy_max = [max(g1.get_xlim()), max(g1.get_ylim())]
                    ax.plot([max(xy_min), min(xy_max)], [max(xy_min), min(xy_max)], c='r')

                # Set axis/title text:
                g1.set_title('Correlation between {} and {}'.format(nam1, nam2))
                # g1.grid(True, axis='y')

                if charge_plot is True:
                    g1.set_xlabel('Charge (%), {}'.format(nam1));
                    g1.set_ylabel('Charge (%), {}'.format(nam2));
                elif log is True:
                    g1.set_xlabel('Abundance (log10 RPM), {}'.format(nam1));
                    g1.set_ylabel('Abundance (log10 RPM), {}'.format(nam2));
                else:
                    g1.set_xlabel('Abundance (RPM), {}'.format(nam1));
                    g1.set_ylabel('Abundance (RPM), {}'.format(nam2));

                # Write to PDF file:
                fig.tight_layout()
                pp.savefig(fig, bbox_inches='tight')
                plt.close(fig)

    def plot_non_temp(self, end='5p', plot_name='_5p-non-template_logo', \
                      n_jobs=4, seq_len_percentile=99, seq_len=None, \
                      _3p_cover=False, verbose=True, plot_log=False):
        if end == '3p':
            col_sele = '3p_non-temp'
        elif end == '5p':
            col_sele = '5p_non-temp'
        else:
            raise Exception('end input has to be either 3p and 5p, not: {}'.format(end))
        self.verbose = verbose
        if verbose:
            print('\nNow collecting data for sample:', end='')

        # Check data exists:
        for _, row in self.sample_df.iterrows():
            stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(stats_fnam))

        # Read each stats CSV file and generate its count matrix:
        data = list()
        for index, row in self.sample_df.copy().iterrows():
            row['end'] = end
            row['col_sele'] = col_sele
            row['seq_len'] = seq_len
            row['seq_len_percentile'] = seq_len_percentile
            row['_3p_cover'] = _3p_cover
            data.append((index, row))
        # Collect data for plotting for each sample:
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self._collect_no_temp_mat, data)

        # Make logo:
        self._run_logomaker(results, plot_name)
        if plot_log:
            # Also make it as log10 transformed:
            self._run_logomaker(results, plot_name, log10=True)

    def _collect_no_temp_mat(self, index, row):
        if self.verbose:
            print('  {}'.format(row['sample_name_unique']), end='')
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)

        col_sele = row['col_sele']
        if row['seq_len'] is None:
            # Limit the length of extracted 5p non-template bases to "seq_len_percentile":
            seq_len = np.percentile([len(s) for s in sample_stats[col_sele]], row['seq_len_percentile'], method='nearest')
        else:
            seq_len = row['seq_len']

        title = 'Sample: {}, logo from {} obs'.format(row['sample_name_unique'], row['N_mapped'])
        mask = np.array([len(s) <= seq_len for s in sample_stats[col_sele].values])
        if row['_3p_cover']:
            mask = mask & (sample_stats['3p_cover'] == True)
        if mask.sum() > 1:
            # The dash "-" sign is not counted in "counts_mat":
            if row['end'] == '3p':
                seq_list = [s + '-'*(seq_len - len(s)) for s in sample_stats.loc[mask, col_sele].values if len(s) <= seq_len]
            else:
                seq_list = ['-'*(seq_len - len(s)) + s for s in sample_stats.loc[mask, col_sele].values if len(s) <= seq_len]
        # Will throw an exception if all gaps:
        try:
            with warnings.catch_warnings(): # ignoring a warning about array assignment
                warnings.simplefilter(action='ignore', category=FutureWarning)
                counts_mat = lm.alignment_to_matrix(seq_list, sample_stats.loc[mask, self.RPM_count_col].values)
        except:
            counts_mat = False
        return([counts_mat, title])

    def plot_UMI_logo(self, plot_name='UMI_logo', n_jobs=4, verbose=True):
        self.verbose = verbose
        if verbose:
            print('\nNow collecting data for sample:', end='')
        # Check files exists before starting:
        for _, row in self.sample_df.iterrows():
            stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(stats_fnam))

        # Read each stats CSV file and generate its count matrix:
        data = list(self.sample_df.iterrows())
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self._collect_UMI_mat, data)
        self._run_logomaker(results, plot_name)

    def _collect_UMI_mat(self, index, row):
        if self.verbose:
            print('  {}'.format(row['sample_name_unique']), end='')
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)
        try:
            UMI_mask = sample_stats['5p_UMI'] != ''
            with warnings.catch_warnings(): # ignoring a warning about array assignment
                warnings.simplefilter(action='ignore', category=FutureWarning)
                counts_mat = lm.alignment_to_matrix(sample_stats[UMI_mask]['5p_UMI'].values)
            title = 'Sample: {}, logo from {} obs (excluding common seqs), {} % of obs vs exp'.format(row['sample_name_unique'], UMI_mask.sum(), round(row['percent_UMI_obs-vs-exp']))
        except:
            counts_mat = False
            title = ''

        return([counts_mat, title])

    def _run_logomaker(self, plot_data, plot_name, log10=False):
        if self.verbose:
            print('\nNow plotting logo plot.', end='')
        if log10:
            plot_name = plot_name + '_log10'
        # Print each to the same PDF file:
        logo_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(logo_fnam) as pp:
            for dat in plot_data:
                counts_mat, title = dat
                if counts_mat is False:
                    continue
                if log10:
                    counts_mat = np.log10(counts_mat+1)
                # logomaker prints a warning when encountering "N"
                # characters, we don't want that:
                with contextlib.redirect_stdout(None):
                    logo_plot = lm.Logo(counts_mat, color_scheme='classic');
                logo_plot.ax.set_title(title, fontsize=15)
                logo_plot.ax.set_xlabel("5' to 3'")
                logo_plot.ax.set_ylabel("Count");
                logo_plot.fig.tight_layout()
                pp.savefig(logo_plot.fig, bbox_inches='tight')
                plt.close(logo_plot.fig)

    def plot_coverage(self, compartment='cyto', plot_type='needle', \
        aa_norm=False, y_norm=False, plot_name='cov_plot_cyto_needle', \
        sample_list=None, verbose=True, n_jobs=4, max_5p_non_temp=10):
        # Check input:
        if compartment != 'mito' and compartment != 'cyto':
            raise Exception('Unknown compartment specified: {}\nValid strings are either either "mito" or "cyto"'.format(compartment))
        if plot_type != 'behrens' and plot_type != 'needle':
            raise Exception('Unknown plot type specified: {}\nValid strings are either either "behrens" or "needle".'.format(plot_type))
        if sample_list is None:
            sample_list = [row['sample_name_unique'] for _, row in self.sample_df.iterrows()]
        sample_list = set(sample_list)

        # Columns used for data aggregation:
        aa_cov_cols = ['sample_name_unique', 'tRNA_annotation_len', \
                       'align_5p_idx', 'align_3p_idx', 'AA_letter', \
                       self.RPM_count_col]
        # Check data exists:
        for _, row in self.sample_df.iterrows():
            stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
            assert(os.path.exists(stats_fnam))
        if verbose:
            print('\nNow collecting data for sample:', end='')

        # Read each stats CSV file and generate its coverage matrix:
        data = list()
        for index, row in self.sample_df.copy().iterrows():
            if not row['sample_name_unique'] in sample_list:
                continue
            row['compartment'] = compartment
            row['aa_norm'] = aa_norm
            row['y_norm'] = y_norm
            data.append((index, row, aa_cov_cols, max_5p_non_temp, verbose))
        # Collect data for plotting for each sample:
        with WorkerPool(n_jobs=n_jobs) as pool:
            results = pool.map(self._collect_coverage_data, data)
        if verbose:
            print('\nNow plotting sample:', end='')

        # Use a 20 color colormap:
        cmap = mpl.colormaps['tab20']
        # Print each plot to the same PDF file:
        cov_fnam = '{}/{}.pdf'.format(self.plotting_dir_abs, plot_name)
        with PdfPages(cov_fnam) as pp:
            # Loop through and generate plots for each sample:
            for res in results:
                if res is False:
                    continue
                cov_count, cov_count_sum, aa_ordered_list, title_info = res
                if verbose:
                    print('  {}'.format(title_info[0]), end='')
                try: # Catch odd errors, like empty dataframe etc. 
                    title = 'Sample: {}, coverage from {} obs'\
                    .format(title_info[0], title_info[1])
                    
                    if plot_type == 'behrens': # Behrens type coverage plot
                        fig, axes = plt.subplots(1, 2, figsize=(13, 7), \
                                                 sharey=False, gridspec_kw={'width_ratios': [1, 7]})

                        # Plot small 10 nt. window from the 5p end:
                        for i in range(1, len(cov_count_sum)):
                            axes[0].stairs(cov_count_sum[i, 0:10], baseline=cov_count_sum[i-1, 0:10], \
                                           fill=True, alpha=0.7, color=cmap(i-1))
                        # Add a black line to outline the top:
                        axes[0].stairs(cov_count_sum[i, 0:10], baseline=None, fill=False, \
                                       edgecolor='black', alpha=1, linewidth=1)

                        # Plot the whole window, from 5p to 3p:
                        for i in range(1, len(cov_count_sum)):
                            axes[1].stairs(cov_count_sum[i], baseline=cov_count_sum[i-1], fill=True, \
                                           alpha=0.7, label=aa_ordered_list[i-1], color=cmap(i-1))
                        axes[1].stairs(cov_count_sum[i], baseline=None, fill=False, \
                                       edgecolor='black', alpha=1, linewidth=1)

                        # Invert the legend labels so they are the same order
                        # as the stacks appear on the plot:
                        handles, labels = axes[1].get_legend_handles_labels()
                        axes[1].legend(handles[::-1], labels[::-1], title='Amino acid', \
                                       bbox_to_anchor=(1.01,1), borderaxespad=0, ncol=2);

                        # Add title and axis labels:
                        if aa_norm:
                            axes[0].set_ylabel("Normalized coverage\n(amino acids equally weighed at 3')")
                        elif y_norm:
                            axes[0].set_ylabel("Normalized coverage (%)")
                        else:
                            axes[0].set_ylabel('Read count')
                        axes[1].set_xlabel("5' to 3' index (mapped to longest tRNA)");
                        axes[1].set_title(title)

                    elif plot_type == 'needle': # Needle-like type coverage plot
                        # Each "needle" follows the midpoint of
                        # the summed count for the 3p (last) position:
                        last_col = cov_count_sum[:, -1]
                        last_col_mid = np.zeros(last_col.shape[0]-1)
                        for i in range(1, last_col.shape[0]):
                            delta = last_col[i] - last_col[i-1]
                            last_col_mid[i-1] = last_col[i-1] + delta/2

                        # Generate top and bottom values at each position
                        # from 5p to 3p for each "needle":
                        cov_funnel_top = np.zeros(cov_count.shape)
                        cov_funnel_bot = np.zeros(cov_count.shape)
                        for i in range(0, cov_count.shape[1]):
                            cov = cov_count[:, i]
                            cov_funnel_top[:, i] = last_col_mid + cov/2
                            cov_funnel_bot[:, i] = last_col_mid - cov/2

                        fig, ax = plt.subplots(1, 1, figsize=(10, 7))
                        for i in range(0, len(cov_funnel_top)):
                            # Plot each "needle":
                            ax.stairs(cov_funnel_top[i, :], baseline=cov_funnel_bot[i, :], fill=True, \
                                      alpha=0.7, label=aa_ordered_list[i], color=cmap(i))

                            # Plot an outline around the needle:
                            funnel_mask = cov_funnel_top[i] > cov_funnel_bot[i]
                            edges = np.where(funnel_mask)[0]
                            edges = np.hstack((edges, edges[-1]+1))
                            ax.stairs(cov_funnel_top[i, funnel_mask], edges, \
                                      baseline=cov_funnel_bot[i, funnel_mask], \
                                      edgecolor='black', alpha=1, linewidth=1)
                        # Invert the legend labels so they are the same order
                        # as the stacks appear on the plot:
                        handles, labels = ax.get_legend_handles_labels()
                        ax.legend(handles[::-1], labels[::-1], title='Amino acid', \
                                  bbox_to_anchor=(1.01,1), borderaxespad=0, ncol=2);

                        # Add title and axis labels:
                        if aa_norm:
                            ax.set_ylabel("Normalized coverage\n(amino acids equally weighed at 3')")
                        elif y_norm:
                            ax.set_ylabel("Normalized coverage (%)")
                        else:
                            ax.set_ylabel('Read count')
                        ax.set_xlabel("5' to 3' index (mapped to longest tRNA)");
                        ax.set_title(title)

                    # Write to PDF file:
                    fig.tight_layout()
                    pp.savefig(fig, bbox_inches='tight')
                    plt.close(fig)
                except Exception as err:
                    pass
                    print(err)

    def _collect_coverage_data(self, index, row, aa_cov_cols, \
                               max_5p_non_temp, verbose):
        print('  {}'.format(row['sample_name_unique']), end='')
        # Read data and add necessary columns:
        stats_fnam = '{}/{}_stats.csv.bz2'.format(self.stats_dir_abs, row['sample_name_unique'])
        with bz2.open(stats_fnam, 'rt', encoding="utf-8") as stats_fh:
            sample_stats = pd.read_csv(stats_fh, keep_default_na=False, dtype=self.stats_csv_header_td)
        # Create single codon filter, to filter out sequences that map to
        # tRNA sequences with different codon/anticodon:
        single_codon = list()
        for anno_str, anticodon in zip(sample_stats['tRNA_annotation'].values, sample_stats['anticodon'].values):
            sc = True
            anno_list = anno_str.split('@')
            for anno in anno_list:
                if not anno.split('-')[2] == anticodon:
                    sc = False
            single_codon.append(sc)
        sample_stats['single_codon'] = single_codon
        sample_stats['Syn_ctr'] = ['Escherichia_coli' in anno and sp != 'ecoli' for sp, anno in zip(sample_stats['species'].values, sample_stats['tRNA_annotation'].values)]
        sample_stats['mito_codon'] = ['mito_tRNA' in anno for anno in sample_stats['tRNA_annotation'].values]
        # Translate codon into single letter amino acid:
        sample_stats['AA_letter'] = [AAA2A[AAA] for AAA in sample_stats['amino_acid'].values]

        # Choose rows from requested compartment:
        type_mask = (sample_stats['3p_cover'] == True) & \
                    (sample_stats['5p_non-temp'].apply(len) <= max_5p_non_temp) & \
                    ((sample_stats['align_3p_nts'] == 'CA') | (sample_stats['align_3p_nts'] == 'CC') | (sample_stats['align_3p_nts'] == 'GA') | (sample_stats['align_3p_nts'] == 'CG')) & \
                    (sample_stats['single_codon']) & \
                    (~sample_stats['Syn_ctr']) & \
                    (sample_stats['AA_letter'].apply(len) == 1) & \
                    (sample_stats['AA_letter'] != 'X')
        if row['compartment'] == 'mito':
            type_mask &= (sample_stats['mito_codon'])
        else:
            type_mask &= (~sample_stats['mito_codon'])

        # Generate dataframe with coverage information.
        # The coverage can be calculated from the count and align_5p_idx.
        cov_df = sample_stats.loc[type_mask, aa_cov_cols].copy()
        cov_df = cov_df.groupby(aa_cov_cols[:-1], as_index=False).agg({self.RPM_count_col: "sum"}).reset_index(drop=True)

        # tRNAs differ in length, so to plot coverage on the same x-axis,
        # we map the coverage of each tRNA to the longest tRNA in the set.
        # Mapping is done using nearest percentile indexing
        # from the queried tRNA to the longest in the set:
        max_len = cov_df['tRNA_annotation_len'].max()
        min_len = cov_df['tRNA_annotation_len'].min()
        if cov_df.empty:
            return(False)
            #return([False, False, False, 'False'])
        len_map_len = dict() # prepare mapping for all tRNAs
        for len_i in range(min_len, max_len+1):
            len_map_len[len_i] = np.percentile(np.arange(max_len), np.linspace(0, 100, len_i), method='nearest')

        # Get the amino acid one letter names sorted and their indexing:
        aa_ordered_list = sorted(set(cov_df['AA_letter']), reverse=True)
        aa_order = {aa: i for i, aa in enumerate(aa_ordered_list)}

        # Add the coverage counts to a matrix:
        cov_count = np.zeros((len(aa_order), max_len))
        # Each row has an amino acid letter, align_5p_idx and a count.
        # Insert these into the matrix:
        for aa_let, _5pi, alen, count in zip(cov_df['AA_letter'], \
                                             cov_df['align_5p_idx'], \
                                             cov_df['tRNA_annotation_len'], \
                                             cov_df[self.RPM_count_col]):
            aa_idx = aa_order[aa_let]
            _5p_idx = _5pi - 1 # shift alignment index to 0 indexing
            _5p_idx_trans = len_map_len[alen][_5p_idx]
            # Amino acids have multiple codons, each with multiple
            # transcripts that can vary in length, thus add instead of assign:
            cov_count[aa_idx, _5p_idx_trans] += count

        # Each count in the matrix repressents a number of reads
        # and their length mapped to the tRNA. Now this is converted
        # to coverage by summation from left to right i.e. 5p to 3p:
        for i in range(cov_count.shape[0]):
            for j in range(cov_count.shape[1]):
                if j > 0:
                    cov_count[i, j] += cov_count[i, j-1]

        # If "aa_norm" is requested the read count for each
        # amino acid is weighted so all amino acids start of with
        # the same coverage at the 3p and normalized so it sums to 1:
        if row['aa_norm']:
            # Weigh each amino acid:
            cov_count = cov_count.T / cov_count[:, -1]
            # Normalize so coverage at 3p is 1:
            cov_count = cov_count.T * (1/cov_count.shape[1])
        # If "y_norm" is requested the y-axis is normalized to 1
        elif row['y_norm']:
            # Normalize so coverage at 3p is summing to 100:
            cov_count = 100 * cov_count / sum(cov_count[:, -1])

        # "cov_count" is summed on each row from left to right,
        # but for plotting we also need to do summation on each
        # column from top to bottom:
        cov_count_sum = copy.deepcopy(cov_count)
        for i in range(cov_count_sum.shape[0]):
            if i > 0:
                cov_count_sum[i] += cov_count_sum[i-1]
        # Add an additional row of zeroes:
        cov_count_sum = np.vstack((np.zeros(cov_count_sum.shape[1]), cov_count_sum))
        title_info = (row['sample_name_unique'], cov_df[self.RPM_count_col].sum())
        return([cov_count, cov_count_sum, aa_ordered_list, title_info])

def freq2ratio(freq):
    '''Convert frequency to ratio.'''
    return(freq / (1 - freq))

def ratio2freq(ratio):
    '''Convert ratio to frequency.'''
    return(ratio / (ratio + 1))



AAA2A = {
'Ala':  'A',
'Arg':  'R',
'Asn':  'N',
'Asp':  'D',
'Cys':  'C',
'Gln':  'Q',
'Glu':  'E',
'Gly':  'G',
'His':  'H',
'Ile':  'I',
'Ile2':  'I',
'Leu':  'L',
'Leu1': 'L',
'Leu2': 'L',
'Lys':  'K',
'Met':  'M',
'Phe':  'F',
'Pro':  'P',
'SeC':  'SeC',
'Ser':  'S',
'Ser1': 'S',
'Ser2': 'S',
'Thr':  'T',
'Trp':  'W',
'Tyr':  'Y',
'Val':  'V',
'eColiLys': 'K',
'eColiThr': 'T',
'iMet': 'M',
'fMet': 'M',
'Sup': 'X',
'Und': 'X'
}



