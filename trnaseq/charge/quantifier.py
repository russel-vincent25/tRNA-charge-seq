"""
tRNA Charge Quantification Module

This module provides standalone functionality for quantifying tRNA charging levels
from alignment statistics CSV files.

Author: Charge-Specialist
"""

import pandas as pd
import numpy as np
import copy
from typing import Optional, Dict, List, Union


class ChargeQuantifier:
    """
    Quantifies tRNA charging levels from alignment statistics.

    This class extracts the charge quantification logic from the existing
    TRNA_plot._get_charge_df() method and provides a standalone interface
    for charge analysis.

    The charging state is determined by the 3' nucleotides:
    - CA: charged (canonical)
    - CC: uncharged (canonical)
    - GA: charged (non-canonical)
    - CG: uncharged (non-canonical)

    Attributes:
        stats_df (pd.DataFrame): Loaded statistics from CSV file
        charge_count_col (str): Column to use for charge calculation ('count' or 'UMIcount')
        RPM_count_col (str): Column to use for RPM calculations ('count' or 'UMIcount')
        excl_align_gap (bool): Exclude alignments with gaps
        excl_09_fmax (bool): Exclude alignments with fmax_score < 0.9
        charge_df (pd.DataFrame): Calculated charge data
        charge_filt (Dict): Filtered charge data by annotation level
    """

    def __init__(self,
                 stats_csv: str,
                 charge_count: str = 'count',
                 RPM_count: str = 'count',
                 excl_align_gap: bool = False,
                 excl_09_fmax: bool = False):
        """
        Initialize the ChargeQuantifier with alignment statistics.

        Args:
            stats_csv: Path to ALL_stats_aggregate.csv file
            charge_count: Column to use for charge calculation ('count' or 'UMIcount')
            RPM_count: Column to use for RPM calculation ('count' or 'UMIcount')
            excl_align_gap: Exclude alignments with gaps
            excl_09_fmax: Exclude alignments with fmax_score < 0.9

        Raises:
            FileNotFoundError: If stats_csv does not exist
            ValueError: If charge_count or RPM_count are invalid
        """
        # Validate count column inputs
        valid_count_cols = ['count', 'UMIcount']
        if charge_count not in valid_count_cols:
            raise ValueError(f'"charge_count" must be one of {valid_count_cols}')
        if RPM_count not in valid_count_cols:
            raise ValueError(f'"RPM_count" must be one of {valid_count_cols}')

        self.charge_count_col = charge_count
        self.RPM_count_col = RPM_count
        self.excl_align_gap = excl_align_gap
        self.excl_09_fmax = excl_09_fmax
        self.charge_filt = dict()

        # Define expected column types for the aggregated stats file
        # Note: The CSV may have either 'align_3p_nt' or 'align_3p_nts'
        self.stats_agg_cols_td = {
            'sample_name_unique': str,
            'sample_name': str,
            'replicate': int,
            'barcode': str,
            'species': str,
            'tRNA_annotation': str,
            'tRNA_annotation_len': int,
            'unique_annotation': bool,
            '5p_cover': bool,
            'align_3p_nts': str,  # Preferred column name (2-char dinucleotide)
            'align_3p_nt': str,   # Legacy column name (single char)
            'codon': str,
            'anticodon': str,
            'amino_acid': str,
            'align_gap': bool,     # Optional: only in full (non-aggregate) stats
            'fmax_score>0.9': bool, # Optional: only in full (non-aggregate) stats
            '3p_cover': bool,      # Optional: only in full (non-aggregate) stats
            'UMIcount': int,       # Optional: only present with UMI data
            'count': int
        }

        # Load the CSV file
        self.stats_df = self._load_stats_csv(stats_csv)

        # Calculate charge data
        self._calculate_charge()

    def _load_stats_csv(self, stats_csv: str) -> pd.DataFrame:
        """
        Load and preprocess the aggregated stats CSV file.

        Args:
            stats_csv: Path to CSV file

        Returns:
            Preprocessed DataFrame

        Raises:
            FileNotFoundError: If file does not exist
        """
        # Read CSV with appropriate dtypes
        df = pd.read_csv(stats_csv, keep_default_na=False)

        # Handle column name discrepancy: align_3p_nt vs align_3p_nts
        if 'align_3p_nt' in df.columns and 'align_3p_nts' not in df.columns:
            df = df.rename(columns={'align_3p_nt': 'align_3p_nts'})

        # Convert dtypes for relevant columns
        for col, dtype in self.stats_agg_cols_td.items():
            if col in df.columns:
                if dtype == bool:
                    df[col] = df[col].astype(bool)
                elif dtype == int:
                    df[col] = df[col].astype(int)
                elif dtype == str:
                    df[col] = df[col].astype(str)

        # Add amino acid single letter code
        df['AA_letter'] = df['amino_acid'].apply(self._amino_acid_to_letter)

        # Add amino acid-codon string
        df['AA_codon'] = df['amino_acid'] + '-' + df['codon']

        # Add tRNA annotation short form
        df['tRNA_anno_short'] = df['tRNA_annotation'].apply(self._shorten_annotation)

        # Create single codon filter
        df['single_codon'] = df.apply(
            lambda row: self._is_single_codon(row['tRNA_annotation'], row['anticodon']),
            axis=1
        )

        # Create single amino acid filter
        df['single_aa'] = df.apply(
            lambda row: self._is_single_aa(row['tRNA_annotation'], row['amino_acid']),
            axis=1
        )

        # Mark mitochondrial and synthetic control tRNAs
        df['mito_codon'] = df['tRNA_annotation'].str.contains('mito_tRNA')
        df['Syn_ctr'] = df.apply(
            lambda row: 'Synthetic' in row['tRNA_annotation'] and row['species'] != 'ecoli',
            axis=1
        )

        return df

    def _amino_acid_to_letter(self, amino_acid: str) -> str:
        """Convert three-letter amino acid code to single letter."""
        # Remove trailing 1/2 from mitochondrial Leu/Ser
        aa = amino_acid[:-1] if amino_acid and amino_acid[-1] in ['1', '2'] else amino_acid

        # Mapping dictionary
        AAA2A = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Ile2': 'I', 'Leu': 'L', 'Leu1': 'L', 'Leu2': 'L', 'Lys': 'K',
            'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'SeC': 'SeC', 'Ser': 'S',
            'Ser1': 'S', 'Ser2': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y',
            'Val': 'V', 'eColiLys': 'K', 'eColiThr': 'T', 'iMet': 'M',
            'fMet': 'M', 'Sup': 'X', 'Und': 'X'
        }

        return AAA2A.get(aa, 'X')

    def _shorten_annotation(self, annotation: str) -> str:
        """Shorten tRNA annotation for display."""
        an_short = list()
        for an_p in annotation.split('@'):
            if '_tRX' in an_p:
                an_short.append('-'.join(an_p.split('-')[1:]) + '-X')
            else:
                an_short.append('-'.join(an_p.split('-')[1:]))
        return '@'.join(an_short)

    def _is_single_codon(self, annotation: str, anticodon: str) -> bool:
        """Check if all annotations have the same anticodon."""
        anno_list = annotation.split('@')
        for anno in anno_list:
            parts = anno.split('-')
            if len(parts) > 2 and parts[2] != anticodon:
                return False
        return True

    def _is_single_aa(self, annotation: str, amino_acid: str) -> bool:
        """Check if all annotations have the same amino acid."""
        anno_list = annotation.split('@')
        for anno in anno_list:
            parts = anno.split('-')
            if len(parts) > 1 and amino_acid not in parts[1]:
                return False
        return True

    def _calculate_charge(self):
        """
        Calculate tRNA charge percentages from alignment statistics.

        This is the core algorithm extracted from TRNA_plot._get_charge_df().
        It processes the align_3p_nts column to determine charging state:
        - CA/GA: charged (amino acid attached)
        - CC/CG: uncharged (no amino acid)

        When fragment classification columns (is_5p_tRF, is_degraded_tRNA) are
        available, they are included in the charge denominator to match the
        TRNA_plot._get_charge_df() formula exactly:
            charge = CA / (CA + CC + 5p_tRF + degraded)
        Otherwise, the simpler formula is used:
            charge = CA / (CA + CC)
        """
        # Check if fragment columns are available
        self._has_fragments = all(
            col in self.stats_df.columns
            for col in ['is_5p_tRF', 'is_degraded_tRNA', 'is_pre_tRNA']
        )

        # Filter rows based on exclusion criteria
        row_mask = (self.stats_df['count'] > 0)  # all True (dummy)
        if self.excl_align_gap:
            if 'align_gap' not in self.stats_df.columns:
                raise ValueError(
                    "excl_align_gap=True requires 'align_gap' column in the data. "
                    f"Available columns: {list(self.stats_df.columns)}")
            row_mask &= (~self.stats_df['align_gap'])
        if self.excl_09_fmax:
            if 'fmax_score>0.9' not in self.stats_df.columns:
                raise ValueError(
                    "excl_09_fmax=True requires 'fmax_score>0.9' column in the data. "
                    f"Available columns: {list(self.stats_df.columns)}")
            row_mask &= (self.stats_df['fmax_score>0.9'])

        charge_df = self.stats_df.loc[row_mask].copy()

        # Determine which count columns are available
        count_cols = [c for c in ['count', 'UMIcount'] if c in charge_df.columns]
        if self.charge_count_col not in charge_df.columns:
            raise ValueError(
                f"charge_count column '{self.charge_count_col}' not found in data. "
                f"Available count columns: {count_cols}")
        if self.RPM_count_col not in charge_df.columns:
            raise ValueError(
                f"RPM_count column '{self.RPM_count_col}' not found in data. "
                f"Available count columns: {count_cols}")

        # Define columns to use for grouping
        stats_agg_cols = [
            'sample_name_unique', 'sample_name', 'replicate', 'barcode', 'species',
            'tRNA_annotation', 'tRNA_anno_short', 'tRNA_annotation_len', 'unique_annotation',
            '5p_cover', '3p_cover', 'align_3p_nts', 'codon', 'anticodon', 'amino_acid',
            'AA_letter', 'AA_codon', 'single_codon', 'single_aa', 'mito_codon', 'Syn_ctr',
        ] + count_cols

        # Filter to only include existing columns
        stats_agg_cols = [col for col in stats_agg_cols if col in charge_df.columns]

        # Rearrange columns so count columns are last
        other_cols = [col for col in stats_agg_cols if col not in count_cols]
        charge_df = charge_df[other_cols + count_cols]

        # Group by all columns except count columns
        agg_dict = {c: "sum" for c in count_cols}
        charge_df = charge_df.groupby(other_cols, as_index=False).agg(
            agg_dict
        ).reset_index(drop=True)

        # Count charge states based on 3' nucleotides
        ct = self.charge_count_col
        nts = charge_df['align_3p_nts']
        charge_df['CA_count'] = np.where(nts == 'CA', charge_df[ct], 0)
        charge_df['CC_count'] = np.where(nts == 'CC', charge_df[ct], 0)
        charge_df['GA_count'] = np.where(nts == 'GA', charge_df[ct], 0)
        charge_df['CG_count'] = np.where(nts == 'CG', charge_df[ct], 0)

        # Derive full-length vs RT drop-off counts from 5p_cover
        self._has_5p_cover = '5p_cover' in charge_df.columns
        if self._has_5p_cover:
            is_full = charge_df['5p_cover'].astype(bool)
            charge_df['full_length_count'] = np.where(is_full, charge_df[ct], 0)
            charge_df['rt_dropoff_count'] = np.where(~is_full, charge_df[ct], 0)

        # Count fragment types when available (matches TRNA_plot._get_charge_df)
        if self._has_fragments:
            charge_df['5p_tRF'] = np.where(
                charge_df.get('is_5p_tRF', False).astype(bool), charge_df[ct], 0)
            charge_df['degraded'] = np.where(
                charge_df.get('is_degraded_tRNA', False).astype(bool), charge_df[ct], 0)
            charge_df['pre_tRNA'] = np.where(
                charge_df.get('is_pre_tRNA', False).astype(bool), charge_df[ct], 0)

        # Group transcripts with different 3p nt (and different 5p_cover)
        # into a single row per tRNA per sample, summing all count columns.
        charge_df_cols = copy.deepcopy(other_cols)
        for drop_col in ['align_3p_nts', '5p_cover']:
            if drop_col in charge_df_cols:
                charge_df_cols.remove(drop_col)
        charge_df_cols.extend(['CA_count', 'CC_count', 'GA_count', 'CG_count'])
        if self._has_5p_cover:
            charge_df_cols.extend(['full_length_count', 'rt_dropoff_count'])
        if self._has_fragments:
            charge_df_cols.extend(['5p_tRF', 'degraded', 'pre_tRNA'])

        # Separate groupby columns from columns to aggregate
        frag_cols = ['5p_tRF', 'degraded', 'pre_tRNA'] if self._has_fragments else []
        cover_cols = ['full_length_count', 'rt_dropoff_count'] if self._has_5p_cover else []
        non_agg_cols = count_cols + ['CA_count', 'CC_count', 'GA_count', 'CG_count'] + cover_cols + frag_cols
        groupby_cols = [col for col in charge_df_cols if col not in non_agg_cols]

        agg_dict2 = {c: "sum" for c in count_cols}
        agg_dict2.update({'CA_count': "sum", 'CC_count': "sum", 'GA_count': "sum", 'CG_count': "sum"})
        if self._has_5p_cover:
            agg_dict2.update({'full_length_count': "sum", 'rt_dropoff_count': "sum"})
        if self._has_fragments:
            agg_dict2.update({'5p_tRF': "sum", 'degraded': "sum", 'pre_tRNA': "sum"})
        charge_df = charge_df.groupby(groupby_cols, as_index=False).agg(
            agg_dict2
        ).reset_index(drop=True)

        # Calculate charge percentages with safe division
        # When fragment data is available, include in denominator (matches TRNA_plot)
        if self._has_fragments:
            ca_cc_frag = charge_df['CA_count'] + charge_df['CC_count'] + charge_df['5p_tRF'] + charge_df['degraded']
            charge_df['charge_canonical'] = np.where(
                ca_cc_frag > 0, 100 * charge_df['CA_count'] / ca_cc_frag, np.nan)
            ga_cg_frag = charge_df['GA_count'] + charge_df['CG_count'] + charge_df['5p_tRF'] + charge_df['degraded']
            charge_df['charge_non-canonical'] = np.where(
                ga_cg_frag > 0, 100 * charge_df['GA_count'] / ga_cg_frag, np.nan)
            charge_df['tRNA_fragment'] = np.where(
                charge_df['count'] > 0,
                (charge_df['5p_tRF'] + charge_df['degraded']) / charge_df['count'], np.nan)
        else:
            ca_cc = charge_df['CA_count'] + charge_df['CC_count']
            charge_df['charge_canonical'] = np.where(
                ca_cc > 0, 100 * charge_df['CA_count'] / ca_cc, np.nan)
            ga_cg = charge_df['GA_count'] + charge_df['CG_count']
            charge_df['charge_non-canonical'] = np.where(
                ga_cg > 0, 100 * charge_df['GA_count'] / ga_cg, np.nan)

        # Calculate RT drop-off fraction when 5p_cover data is available
        if self._has_5p_cover:
            total = charge_df['full_length_count'] + charge_df['rt_dropoff_count']
            charge_df['rt_dropoff_fraction'] = np.where(
                total > 0,
                charge_df['rt_dropoff_count'] / total,
                np.nan
            )

        # Add sample total count for RPM calculation
        df_count = charge_df[~charge_df['Syn_ctr']].groupby(
            ['sample_name_unique'],
            as_index=False
        ).agg({self.RPM_count_col: "sum"}).reset_index(drop=True)

        charge_df = charge_df.merge(df_count, on='sample_name_unique', suffixes=('', '_sample_tot'))

        # Calculate RPM
        charge_df['RPM'] = charge_df[self.RPM_count_col] / (
            charge_df[f'{self.RPM_count_col}_sample_tot'] / 1e6
        )
        charge_df = charge_df.drop(columns=[f'{self.RPM_count_col}_sample_tot'])

        # Store the main charge dataframe
        self.charge_df = charge_df

        # Create filtered versions by annotation level
        self._create_filtered_charge_dfs()

    def _create_filtered_charge_dfs(self):
        """Create filtered charge dataframes by amino acid, codon, and transcript."""
        # Build dynamic agg dict based on available count columns
        count_cols = [c for c in ['count', 'UMIcount'] if c in self.charge_df.columns]
        base_agg = {c: "sum" for c in count_cols}
        base_agg.update({
            'CA_count': "sum", 'CC_count': "sum",
            'GA_count': "sum", 'CG_count': "sum",
            'RPM': "sum"
        })
        if self._has_5p_cover:
            base_agg.update({'full_length_count': "sum", 'rt_dropoff_count': "sum"})
        if self._has_fragments:
            base_agg.update({'5p_tRF': "sum", 'degraded': "sum", 'pre_tRNA': "sum"})

        def _recalc_charge(df):
            """Recalculate charge percentages and RT drop-off after aggregation."""
            if self._has_fragments:
                ca_cc_frag = df['CA_count'] + df['CC_count'] + df['5p_tRF'] + df['degraded']
                df['charge_canonical'] = np.where(
                    ca_cc_frag > 0, 100 * df['CA_count'] / ca_cc_frag, np.nan)
                ga_cg_frag = df['GA_count'] + df['CG_count'] + df['5p_tRF'] + df['degraded']
                df['charge_non-canonical'] = np.where(
                    ga_cg_frag > 0, 100 * df['GA_count'] / ga_cg_frag, np.nan)
                df['tRNA_fragment'] = np.where(
                    df['count'] > 0, (df['5p_tRF'] + df['degraded']) / df['count'], np.nan)
            else:
                ca_cc = df['CA_count'] + df['CC_count']
                df['charge_canonical'] = np.where(
                    ca_cc > 0, 100 * df['CA_count'] / ca_cc, np.nan)
                ga_cg = df['GA_count'] + df['CG_count']
                df['charge_non-canonical'] = np.where(
                    ga_cg > 0, 100 * df['GA_count'] / ga_cg, np.nan)
            if self._has_5p_cover:
                total = df['full_length_count'] + df['rt_dropoff_count']
                df['rt_dropoff_fraction'] = np.where(
                    total > 0, df['rt_dropoff_count'] / total, np.nan)
            return df

        # Filter by amino acid
        aa_mask = self.charge_df['single_aa']
        charge_df_aa = self.charge_df[aa_mask].groupby([
            'sample_name_unique', 'sample_name', 'replicate', 'barcode',
            'amino_acid', 'AA_letter', 'mito_codon', 'Syn_ctr'
        ], as_index=False).agg(base_agg).reset_index(drop=True)
        self.charge_filt['aa'] = _recalc_charge(charge_df_aa)

        # Filter by codon
        cd_mask = self.charge_df['single_codon']
        charge_df_cd = self.charge_df[cd_mask].groupby([
            'sample_name_unique', 'sample_name', 'replicate', 'barcode',
            'codon', 'anticodon', 'AA_codon', 'amino_acid', 'AA_letter',
            'mito_codon', 'Syn_ctr'
        ], as_index=False).agg(base_agg).reset_index(drop=True)
        self.charge_filt['codon'] = _recalc_charge(charge_df_cd)

        # Filter by transcript
        tr_mask = self.charge_df['unique_annotation']
        charge_df_tr = self.charge_df[tr_mask].groupby([
            'sample_name_unique', 'sample_name', 'replicate', 'barcode',
            'tRNA_annotation', 'tRNA_anno_short', 'tRNA_annotation_len',
            'codon', 'anticodon', 'AA_codon', 'amino_acid', 'AA_letter',
            'mito_codon', 'Syn_ctr'
        ], as_index=False).agg(base_agg).reset_index(drop=True)
        self.charge_filt['tr'] = _recalc_charge(charge_df_tr)

    def quantify_all(self,
                     level: str = 'transcript',
                     include_synthetic: bool = False,
                     include_mito: bool = True) -> pd.DataFrame:
        """
        Calculate charge percentages for all tRNAs across all samples.

        Args:
            level: Annotation level ('transcript', 'codon', or 'aa')
            include_synthetic: Include synthetic control tRNAs
            include_mito: Include mitochondrial tRNAs

        Returns:
            DataFrame with charge data including:
            - sample_name_unique, sample_name, replicate
            - tRNA identifiers (depends on level)
            - CA_count, CC_count, GA_count, CG_count
            - charge_canonical, charge_non-canonical
            - RPM, count, UMIcount

        Raises:
            ValueError: If level is not valid
        """
        if level == 'aa':
            df = self.charge_filt['aa'].copy()
        elif level == 'codon':
            df = self.charge_filt['codon'].copy()
        elif level == 'transcript':
            df = self.charge_filt['tr'].copy()
        else:
            raise ValueError(f'Unknown level: {level}. Must be "aa", "codon", or "transcript".')

        # Apply filters
        if not include_synthetic:
            df = df[~df['Syn_ctr']]
        if not include_mito:
            df = df[~df['mito_codon']]

        return df.reset_index(drop=True)

    def quantify_single_trna(self,
                            trna_id: str,
                            level: str = 'transcript') -> pd.DataFrame:
        """
        Calculate charge for a specific tRNA across all samples.

        Args:
            trna_id: tRNA identifier (depends on level):
                - transcript: full annotation (e.g., 'Homo_sapiens_tRNA-Ala-AGC-1-1')
                - codon: AA_codon format (e.g., 'Ala-GCT')
                - aa: amino acid name (e.g., 'Ala')
            level: Annotation level ('transcript', 'codon', or 'aa')

        Returns:
            DataFrame with charge data for the specified tRNA

        Raises:
            ValueError: If level is not valid or tRNA not found
        """
        df = self.quantify_all(level=level, include_synthetic=True, include_mito=True)

        # Filter by identifier based on level
        if level == 'transcript':
            result = df[df['tRNA_annotation'].str.contains(trna_id, case=False, regex=False)]
        elif level == 'codon':
            result = df[df['AA_codon'] == trna_id]
        elif level == 'aa':
            result = df[df['amino_acid'] == trna_id]
        else:
            raise ValueError(f'Unknown level: {level}')

        if result.empty:
            raise ValueError(f'tRNA "{trna_id}" not found at level "{level}"')

        return result.reset_index(drop=True)

    def get_summary_statistics(self, level: str = 'transcript') -> pd.DataFrame:
        """
        Get summary statistics of charge across samples.

        Args:
            level: Annotation level ('transcript', 'codon', or 'aa')

        Returns:
            DataFrame with mean, std, min, max charge by tRNA
        """
        df = self.quantify_all(level=level)

        # Group by tRNA identifier
        if level == 'transcript':
            group_col = 'tRNA_annotation'
        elif level == 'codon':
            group_col = 'AA_codon'
        else:
            group_col = 'amino_acid'

        summary = df.groupby(group_col).agg({
            'charge_canonical': ['mean', 'std', 'min', 'max', 'count'],
            'RPM': ['mean', 'std', 'min', 'max']
        }).reset_index()

        # Flatten column names
        summary.columns = [
            '_'.join(col).strip('_') if col[1] else col[0]
            for col in summary.columns.values
        ]

        return summary

    def export_to_csv(self,
                     output_file: str,
                     level: str = 'transcript',
                     include_synthetic: bool = False,
                     include_mito: bool = True):
        """
        Export charge data to CSV file.

        Args:
            output_file: Path to output CSV file
            level: Annotation level ('transcript', 'codon', or 'aa')
            include_synthetic: Include synthetic control tRNAs
            include_mito: Include mitochondrial tRNAs
        """
        df = self.quantify_all(
            level=level,
            include_synthetic=include_synthetic,
            include_mito=include_mito
        )
        df.to_csv(output_file, index=False)
        print(f'Exported charge data to {output_file}')
