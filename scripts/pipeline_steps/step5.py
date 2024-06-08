import pandas as pd
from scripts.utils.sequence_utils import get_nucleotide_at_position
from scripts.pyensembl_operations import *
from scripts.globals import *
import os
import sqlite3
from scripts.globals import MUTSIG_PROBABILITIES


def split_id_column(df):
    df[['vcf_id', 'chr', 'pos', 'ref', 'alt', 'mirna_accession']
       ] = df['id'].str.split('_', expand=True)
    df.pos = df.pos.astype(int)
    return df


def create_mutation_context_string(row):
    try:
        ref = row['ref']
        alt = row['alt']
        before = row['before']
        after = row['after']

        if ref in ['C', 'T']:  # Pyrimidine
            return f"{before}[{ref}>{alt}]{after}"
        elif ref in ['A', 'G']:  # Purine
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            ref_complement = complement[ref]
            alt_complement = complement[alt]
            before_complement = ''.join(
                complement[base] for base in before[::-1])
            after_complement = ''.join(
                complement[base] for base in after[::-1])
            return f"{after_complement}[{ref_complement}>{alt_complement}]{before_complement}"
        else:
            raise ValueError(f"Invalid reference nucleotide: {ref}")
    except KeyError as e:
        print(f"Error: Missing column - {str(e)}")
        return ''
    except Exception as e:
        print(f"Error: {str(e)}")
        return ''


def generate_mutation_context_column(df):
    df["before"] = df.apply(
        lambda x: get_nucleotide_at_position(x['chr'], x["pos"]-1), axis=1)
    df["after"] = df.apply(lambda x: get_nucleotide_at_position(
        x['chr'], x["pos"]+1), axis=1)
    df['mutation_context'] = df.apply(create_mutation_context_string, axis=1)
    df["mutsig_key"] = df["vcf_id"] + "_" + df["mutation_context"]
    df.drop(columns=["before", "after", "ref", "alt"], inplace=True)
    return df


def add_mutsig_probabilities(df, mutsig_file):
    mutsig_df = pd.read_csv(mutsig_file)
    mutsig_df["mutsig_key"] = mutsig_df["Sample Names"] + \
        "_" + mutsig_df["MutationTypes"]
    mutsig_dict = mutsig_df.set_index('mutsig_key')['mutsig'].to_dict()
    df["mutsig"] = df["mutsig_key"].map(mutsig_dict)
    return df



def generate_is_intron_column(df, assembly):

    @lru_cache(maxsize=None)
    def cached_pyensembl_call(locus):

        chrom, pos = locus.split(':')
        result = assembly.exons_at_locus(chrom, int(pos))

        return np.nan if len(result) == 0 else result[0].exon_id

    unique_loci = df['locus'].unique()
    map_dict = {locus: cached_pyensembl_call(locus) for locus in unique_loci}
    df['is_exon'] = df['locus'].map(map_dict)
    
    mask = (df.is_exon.isna()) & (~df.gene_id.isna())

    df["is_intron"] = mask

    df.drop(columns=["is_exon"], inplace=True)

    return df



###################################################


def generate_gene_id_column(df, assembly):
    unique_loci = df['locus'].unique()
    
    @lru_cache(maxsize=None)
    def cached_pyensembl_call(locus):

        chrom, pos = locus.split(':')
        result = assembly.gene_ids_at_locus(chrom, int(pos))

        return np.nan if len(result) == 0 else result[0]

    gene_id_dict = {locus: cached_pyensembl_call(locus) for locus in unique_loci}
    df['gene_id'] = df['locus'].map(gene_id_dict)

    return df


def filter_rows_with_same_prediction(df, threshold=0.5):
    """
    Filter rows from a DataFrame based on the difference between the wt and mut predictions.

    Args:
        df (pandas.DataFrame): The input DataFrame.
        threshold (float, optional): The threshold to determine if the difference is significant. Default is 0.5.

    Returns:
        pandas.DataFrame: The filtered DataFrame.
    """
    # Create a mask based on the difference between wt and mut predictions.
    # The mask is True if the difference is significant (i.e., the predictions are different), and False otherwise.
    mask = (df['wt_prediction'] >= threshold) != (df['mut_prediction'] >= threshold)
    
    # Filter the DataFrame based on the mask and reset the index.
    df = df[mask]
    
    # Reset the index and return the filtered DataFrame.
    return df.reset_index(drop=True)

    


def apply_step_5(file_path, assembly, mutsig_probabilities):
    
    df = pd.read_csv(file_path)
    df = filter_rows_with_same_prediction(df)
    
    df = split_id_column(df)
    df["locus"] = df["chr"] + ":" + df["pos"].astype(str)
    
    df = generate_gene_id_column(df, assembly)
    df = generate_mutation_context_column(df)
    df = add_mutsig_probabilities(df, mutsig_probabilities)
    df = generate_is_intron_column(df, assembly)
    
    df["is_gain"] = df.mut_prediction > df.wt_prediction
    df["is_gene_upregulated"] = ~df.is_gain

    
    df.drop(columns=["chr", "pos", "locus", "mutsig_key"], inplace=True)

    return df

def crawl_and_import_results(folder_path, ending_string, db_path, table_name, assembly):
    csv_files = []

    # Find CSV files
    for root, _, files in os.walk(folder_path):
        csv_files.extend(
            os.path.join(root, file)
            for file in files
            if file.endswith(f"{ending_string}.csv")
        )

    # Connect to SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create table if it doesn't exist
    cursor.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            id TEXT PRIMARY KEY,
            wt_prediction REAL,
            mut_prediction REAL,
            pred_difference REAL,
            vcf_id TEXT,
            mirna_accession TEXT,
            gene_id TEXT,
            mutation_context TEXT,
            mutsig TEXT,
            is_intron BOOLEAN,
            is_gain BOOLEAN,
            is_gene_upregulated BOOLEAN
        )
    """)

    # Import CSV files into the table
    for file_path in csv_files:
        
        df = apply_step_5(file_path, assembly, MUTSIG_PROBABILITIES)

        columns = ', '.join(df.columns)
        placeholders = ', '.join(['?'] * len(df.columns))

        # Insert data into the table
        for row in df.values:
            cursor.execute(f"INSERT INTO {table_name} ({columns}) VALUES ({placeholders})", row)

    # Commit changes and close the connection
    conn.commit()
    conn.close()
