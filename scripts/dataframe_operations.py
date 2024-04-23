import pandas as pd
from scripts.globals import TA_SPS_CSV
from scripts.sequence_operations import calculate_au_content

def create_results_dataframe(wt_array, mutated_array):
    """Create a consolidated dataframe from the results."""
    wt_result_df = pd.DataFrame(wt_array)
    mut_result_df = pd.DataFrame(mutated_array)
    df = pd.concat([wt_result_df, mut_result_df])
    
    colnames = ["mrna_start", "mrna_end", "mrna_dot_bracket_5to3", "mirna_start", "mirna_end", "mirna_dot_bracket_5to3", "pred_energy", "mutation_id", "mirna_accession", "mrna_sequence", "mirna_sequence", "is_mutated"]
    df.columns = colnames
    
    df["id"] = df["mutation_id"].astype(str) + "_" + df["mirna_accession"].astype(str)
    df.drop(columns=["mutation_id"], inplace=True)
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
    mid_strings = df['mirna_dot_bracket_5to3'].apply(lambda x: ''.join('1' if char == ')' else '0' for char in x))
    end_strings = (df['mirna_sequence'].str.len() - df['mirna_end'] - 1).apply(lambda x: '0' * x)

    df['alignment_string'] = start_strings + mid_strings + end_strings

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
    df["pred_num_basepairs"] = df["alignment_string"].str.count("1")

    # Count the number of '1' characters in the seed region (positions 2-7)
    df["pred_seed_basepairs"] = df["alignment_string"].str.slice(1, 7).str.count("1")

    return df



def generate_ta_sps_columns(df):
    """
    Add 'ta_log10' and 'sps_mean' columns to the input DataFrame based on the miRNA seed sequence.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mirna_sequence' column.

    Returns:
        pandas.DataFrame: The input DataFrame with 'ta_log10' and 'sps_mean' columns added.
    """
    # Generate temporary seed column
    df["seed"] = df["mirna_sequence"].str[1:8].str.replace("T", "U")
    # Read ta sps data
    ta_sps_df = pd.read_csv(TA_SPS_CSV, usecols=["seed_8mer", "ta_log10", "sps_mean"])
    ta_sps_df = ta_sps_df.rename(columns={"seed_8mer": "seed"})
    # Merge dataframes on seed column
    df = df.merge(ta_sps_df, on="seed", how="left")
    # Drop temporary column
    df.drop(columns=["seed"], inplace=True)

    return df


def generate_important_sites(df):
    """
    Generate columns for important miRNA target site features.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mre_region' and 'alignment_string' columns.

    Returns:
        pandas.DataFrame: The input DataFrame with additional columns for important site features.
    """
    # Anchor site
    df["anchor_a"] = (df["mre_region"].str[-1] == "A").astype(int)

    # Seed match features
    df["6mer_seed"] = (df["alignment_string"].str[1:7].str.count("0") == 0).astype(int)
    df["match_8"] = (df["alignment_string"].str[7] == "1").astype(int)
    df["6mer_seed_1_mismatch"] = (df["alignment_string"].str[1:7].str.count("0") == 1).astype(int)
    df["empty_seed"] = (df["alignment_string"].str[1:8].str.count("1") == 0).astype(int)

    # Compensatory and supplementary sites
    df["compensatory_site"] = (df["alignment_string"].str[12:17].str.count("0") == 0).astype(int)
    df["supplementary_site"] = (df["alignment_string"].str[12:16].str.count("0") == 0).astype(int)
    df["supplementary_site_2"] = (df["alignment_string"].str[16:21].str.count("0") == 0).astype(int)

    # Consecutive match
    df["9_consecutive_match_anywhere"] = (df["alignment_string"].str.contains("1{" + str(9) + ",}")).astype(int)

    return df


def generate_seed_type_columns(df):
    """
    Generate columns for different types of miRNA seed matches.

    Args:
        df (pandas.DataFrame): A DataFrame containing columns for various miRNA target site features.

    Returns:
        pandas.DataFrame: The input DataFrame with additional columns for seed match types.
    """
    # Canonical seed matches
    has_anchor_a = (df['anchor_a'] == 1)
    has_6mer_seed = (df['6mer_seed'] == 1)
    has_match_8 = (df['match_8'] == 1)
    no_supplementary_sites = (df['supplementary_site'] == 0) & (df['supplementary_site_2'] == 0)

    df['seed_8mer'] = (has_anchor_a & has_6mer_seed & has_match_8).astype(int)
    df['seed_7mer_a1'] = (has_anchor_a & has_6mer_seed & ~has_match_8).astype(int)
    df['seed_7mer_m8'] = (~has_anchor_a & has_6mer_seed & has_match_8 & no_supplementary_sites).astype(int)

    # Non-canonical seed matches
    has_compensatory_site = (df['compensatory_site'] == 1)
    has_6mer_seed_1_mismatch = (df['6mer_seed_1_mismatch'] == 1)
    has_supplementary_site = (df['supplementary_site'] == 1)
    has_supplementary_site_2 = (df['supplementary_site_2'] == 1)
    has_empty_seed = (df['empty_seed'] == 1)
    has_9_consecutive_match = (df['9_consecutive_match_anywhere'] == 1)
    has_many_basepairs = (df['pred_num_basepairs'] > 10)

    df['seed_compensatory'] = (has_compensatory_site & has_6mer_seed_1_mismatch & has_match_8).astype(int)
    df['seed_clash_2'] = (has_supplementary_site & has_6mer_seed & has_match_8).astype(int)
    df['seed_clash_3'] = (has_supplementary_site_2 & has_6mer_seed & has_match_8).astype(int)
    df['seed_clash_4'] = (has_empty_seed & has_9_consecutive_match).astype(int)
    df['seed_clash_5'] = (has_many_basepairs & ~has_6mer_seed).astype(int)

    return df


def generate_local_au_content_column(vcf_df):

    vcf_df["local_au_content"] = vcf_df['mrna_sequence'].apply(calculate_au_content)
    
    return vcf_df