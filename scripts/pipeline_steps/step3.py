
# df = process_rnaduplex_output(rnaduplex_output_file)

# # Step 3: Prediction
# df = process_df_for_prediction(df)
# df, id_array = reorder_columns_for_prediction(df)
# predictions = make_predictions_with_xgb(df)
# df = create_results_df(id_array, predictions, filter_range=QUANTILE_RANGE)

import numpy as np
import pandas as pd
from scripts.globals import *
from functools import lru_cache
import re


def process_rnaduplex_output(rnaduplex_output_file):
    colnames = ["mutation_id", "mirna_accession", "mrna_dot_bracket_5to3", "mirna_dot_bracket_5to3",
                "mrna_start", "mrna_end", "mirna_start", "mirna_end", "pred_energy", "is_mutated"]
    df = pd.read_csv(rnaduplex_output_file, header=None, names=colnames)
    df["id"] = df["mutation_id"] + "_" + \
        df["mirna_accession"] + "_" + df["is_mutated"]
    df = df.sort_values(by=['id', 'is_mutated'], ascending=[False, True])

    # Add miRNA sequence
    mirna_dict = pd.read_csv(MIRNA_CSV).set_index(
        'mirna_accession')['sequence'].to_dict()
    df["mirna_sequence"] = df["mirna_accession"].map(mirna_dict)
    return df


def generate_mirna_conservation_column(df):
    """
    Add a 'mirna_conservation' column to the input DataFrame based on miRNA conservation data.
    Automatically downcast the 'mirna_conservation' column to the most appropriate numerical dtype.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mirna_accession' column.

    Returns:
        pandas.DataFrame: The input DataFrame with a 'mirna_conservation' column added and automatically downcasted.
    """
    mirna_df = (pd.read_csv(MIRNA_CSV,
                            usecols=["mirna_accession", "conservation"])
                .rename(columns={"conservation": "mirna_conservation"})
                [["mirna_accession", "mirna_conservation"]])

    df = df.merge(mirna_df, on="mirna_accession", how="left")

    # Downcast the 'mirna_conservation' column to 'integer' since all values are integers
    df['mirna_conservation'] = pd.to_numeric(
        df['mirna_conservation'], downcast='integer')

    return df


@lru_cache(maxsize=None)
def split_mutation_id_cached(mutation_id):
    return mutation_id.split('_')


def split_mutation_ids(df):
    # Create a list of column names for the new columns
    new_column_names = ['vcf_id', 'chr', 'pos', 'ref', 'alt']

    # Add the new columns to the original DataFrame
    df[new_column_names] = pd.DataFrame(df['mutation_id'].apply(
        split_mutation_id_cached).tolist(), index=df.index, columns=new_column_names)

    # Drop the 'mutation_id' column
    df.drop(['mutation_id', "vcf_id"], axis=1, inplace=True)

    return df


def generate_mre_sequence_column(df):
    """
    Generate the miRNA response element (MRE) sequence for each row in the input DataFrame.

    Args:
        df (pandas.DataFrame): A DataFrame containing columns 'mrna_sequence', 'mrna_end', 'mirna_start', and 'mirna_sequence'.

    Returns:
        pandas.DataFrame: The input DataFrame with new columns 'mre_start', 'mre_end', and 'mre_region' added.
    """
    # Calculate miRNA length
    df["mirna_length"] = df["mirna_sequence"].str.len()

    # Calculate MRE coordinates
    df["mre_end"] = df["mrna_end"] + df["mirna_start"]
    df["mre_start"] = df["mre_end"] - df["mirna_length"]

    # Ensure MRE start is not negative
    df["mre_start"] = df["mre_start"].clip(lower=0)

    # Extract MRE sequence using list comprehension for better performance
    df["mre_region"] = [mrna_seq[start:end] for mrna_seq, start,
                        end in zip(df["mrna_sequence"], df["mre_start"], df["mre_end"])]

    # Drop temporary column
    df.drop(columns=["mirna_length"], inplace=True)

    return df


@lru_cache(maxsize=None)
def calculate_au_content(sequence):
    au_count = sequence.count(
        'A') + sequence.count('T') + sequence.count('U')
    return None if len(sequence) == 0 else au_count / len(sequence)


def generate_local_au_content_column(df):
    # Apply the cached calculate_au_content function to each mrna_sequence in the DataFrame
    df["local_au_content"] = df['mrna_sequence'].apply(calculate_au_content)
    return df


def generate_ta_sps_columns(df):
    """
    Add 'ta_log10' and 'sps_mean' columns to the input DataFrame based on the miRNA seed sequence.
    Downcast the new columns to the smallest numerical dtype possible.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mirna_sequence' column.

    Returns:
        pandas.DataFrame: The input DataFrame with 'ta_log10' and 'sps_mean' columns added and downcasted.
    """
    # Generate temporary seed column
    df["seed"] = df["mirna_sequence"].str.slice(
        1, 8).replace({'T': 'U'}, regex=True)
    # Read ta sps data
    ta_sps_df = pd.read_csv(TA_SPS_CSV, usecols=[
                            "seed_8mer", "ta_log10", "sps_mean"])
    ta_sps_df = ta_sps_df.rename(columns={"seed_8mer": "seed"})
    # Merge dataframes on seed column
    df = df.merge(ta_sps_df, on="seed", how="left")
    # Downcast the new columns
    df['ta_log10'] = pd.to_numeric(df['ta_log10'], downcast='float')
    df['sps_mean'] = pd.to_numeric(df['sps_mean'], downcast='float')
    # Drop temporary column
    df.drop(columns=["seed"], inplace=True)

    return df


def generate_alignment_string_from_dot_bracket(df):
    """
    Generate an alignment string for each row in the DataFrame based on the miRNA dot-bracket structure.

    Args:
        df (pandas.DataFrame): A DataFrame containing columns 'mirna_start', 'mirna_dot_bracket_5to3', 'mirna_sequence', and 'mirna_end'.

    Returns:
        pandas.DataFrame: The input DataFrame with a new column 'alignment_string' added.
    """
    start_strings = df['mirna_start'].apply(lambda x: '0' * x)
    mid_strings = df['mirna_dot_bracket_5to3'].apply(
        lambda x: ''.join('1' if char == ')' else '0' for char in x))
    end_strings = (df['mirna_sequence'].str.len() -
                   df['mirna_end'] - 1).apply(lambda x: '0' * x)

    df['alignment_string'] = start_strings + mid_strings + end_strings
    df.drop(columns=["mrna_dot_bracket_5to3",
            "mirna_dot_bracket_5to3"], inplace=True)

    return df


def generate_match_count_columns(df):
    """
    Generate two new columns in the DataFrame: 'pred_num_basepairs' and 'pred_seed_basepairs'.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'alignment_string' column.

    Returns:
        pandas.DataFrame: The input DataFrame with two new columns added.
    """
    # Count the number of '1' characters in the entire alignment string
    df["pred_num_basepairs"] = df["alignment_string"].str.count(
        "1").astype("uint8")

    # Count the number of '1' characters in the seed region (positions 2-7)
    df["pred_seed_basepairs"] = df["alignment_string"].str.slice(
        1, 7).str.count("1").astype("uint8")

    return df


def generate_important_sites_column(df):
    """
    Optimized function to generate columns for important miRNA target site features.
    Downcast binary columns to the smallest integer dtype.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mre_region' and 'alignment_string' columns.

    Returns:
        pandas.DataFrame: The input DataFrame with additional columns for important site features.
    """
    # Precompile the regular expression for consecutive matches
    consecutive_match_re = re.compile("1{9,}")

    # Anchor site
    df["anchor_a"] = df["mre_region"].str.endswith("A").astype(np.int8)

    # Pre-calculate slices of 'alignment_string' to avoid repeated operations
    alignment_slice_1_7 = df["alignment_string"].str[1:7]
    alignment_slice_7 = df["alignment_string"].str[7]
    alignment_slice_1_8 = df["alignment_string"].str[1:8]
    alignment_slice_12_17 = df["alignment_string"].str[12:17]
    alignment_slice_12_16 = df["alignment_string"].str[12:16]
    alignment_slice_16_21 = df["alignment_string"].str[16:21]

    # Seed match features
    df["6mer_seed"] = alignment_slice_1_7.apply(
        lambda x: x.count("0") == 0).astype(np.int8)
    df["match_8"] = (alignment_slice_7 == "1").astype(np.int8)
    df["6mer_seed_1_mismatch"] = alignment_slice_1_7.apply(
        lambda x: x.count("0") == 1).astype(np.int8)
    df["empty_seed"] = alignment_slice_1_8.apply(
        lambda x: x.count("1") == 0).astype(np.int8)

    # Compensatory and supplementary sites
    df["compensatory_site"] = alignment_slice_12_17.apply(
        lambda x: x.count("0") == 0).astype(np.int8)
    df["supplementary_site"] = alignment_slice_12_16.apply(
        lambda x: x.count("0") == 0).astype(np.int8)
    df["supplementary_site_2"] = alignment_slice_16_21.apply(
        lambda x: x.count("0") == 0).astype(np.int8)

    # Consecutive match
    df["9_consecutive_match_anywhere"] = df["alignment_string"].apply(
        lambda x: bool(consecutive_match_re.search(x))).astype(np.int8)

    return df


def generate_seed_type_columns(df):
    """
    Generate columns for different types of miRNA seed matches.
    Downcast new binary columns to the smallest integer dtype.

    Args:
        df (pandas.DataFrame): A DataFrame containing columns for various miRNA target site features.

    Returns:
        pandas.DataFrame: The input DataFrame with additional columns for seed match types.
    """
    # Canonical seed matches
    has_anchor_a = (df['anchor_a'] == 1)
    has_6mer_seed = (df['6mer_seed'] == 1)
    has_match_8 = (df['match_8'] == 1)
    no_supplementary_sites = (df['supplementary_site'] == 0) & (
        df['supplementary_site_2'] == 0)

    df['seed_8mer'] = (has_anchor_a & has_6mer_seed &
                       has_match_8).astype(np.int8)
    df['seed_7mer_a1'] = (has_anchor_a & has_6mer_seed &
                          ~has_match_8).astype(np.int8)
    df['seed_7mer_m8'] = (~has_anchor_a & has_6mer_seed &
                          has_match_8 & no_supplementary_sites).astype(np.int8)

    # Non-canonical seed matches
    has_compensatory_site = (df['compensatory_site'] == 1)
    has_6mer_seed_1_mismatch = (df['6mer_seed_1_mismatch'] == 1)
    has_supplementary_site = (df['supplementary_site'] == 1)
    has_supplementary_site_2 = (df['supplementary_site_2'] == 1)
    has_empty_seed = (df['empty_seed'] == 1)
    has_9_consecutive_match = (df['9_consecutive_match_anywhere'] == 1)
    has_many_basepairs = (df['pred_num_basepairs'] > 10)

    df['seed_compensatory'] = (
        has_compensatory_site & has_6mer_seed_1_mismatch & has_match_8).astype(np.int8)
    df['seed_clash_2'] = (has_supplementary_site &
                          has_6mer_seed & has_match_8).astype(np.int8)
    df['seed_clash_3'] = (has_supplementary_site_2 &
                          has_6mer_seed & has_match_8).astype(np.int8)
    df['seed_clash_4'] = (
        has_empty_seed & has_9_consecutive_match).astype(np.int8)
    df['seed_clash_5'] = (has_many_basepairs & ~has_6mer_seed).astype(np.int8)

    return df


def generate_mre_au_content_column(df):
    # Apply the cached calculate_au_content function to each mre_region in the DataFrame
    df["mre_au_content"] = df['mre_region'].apply(calculate_au_content)
    return df
