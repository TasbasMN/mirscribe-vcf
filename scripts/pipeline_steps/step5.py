import pandas as pd
from scripts.utils.sequence_utils import get_nucleotide_at_position
from scripts.pyensembl_operations import *
from scripts.globals import *
import os
import sqlite3
from scripts.globals import MUTSIG_PROBABILITIES


def filter_rows_with_same_prediction(df, threshold=0.5):
    mask = (df['wt_prediction'] >= threshold) != (df['mut_prediction'] >= threshold)
    return df[mask].reset_index(drop=True)



def split_id_column(df):
    df[['vcf_id', 'chr', 'pos', 'ref', 'alt', 'mirna_accession']] = df['id'].str.split('_', expand=True)
    df['pos'] = pd.to_numeric(df['pos'], downcast='integer')
    return df



COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def create_mutation_context_string(row):
    ref, alt, before, after = row['ref'], row['alt'], row['before'], row['after']
    if ref in 'CT':  # Pyrimidine
        return f"{before}[{ref}>{alt}]{after}"
    elif ref in 'AG':  # Purine
        ref_complement = COMPLEMENT[ref]
        alt_complement = COMPLEMENT[alt]
        before_complement = ''.join(COMPLEMENT[base] for base in before[::-1])
        after_complement = ''.join(COMPLEMENT[base] for base in after[::-1])
        return f"{after_complement}[{ref_complement}>{alt_complement}]{before_complement}"
    else:
        raise ValueError(f"Invalid reference nucleotide: {ref}")



def generate_mutation_context_column(df):
    df['before'] = df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']-1), axis=1)
    df['after'] = df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']+1), axis=1)
    df['mutation_context'] = df.apply(create_mutation_context_string, axis=1)
    df['mutsig_key'] = df['vcf_id'].astype(str) + '_' + df['mutation_context']
    df.drop(columns=['before', 'after', 'ref', 'alt'], inplace=True)
    return df



def add_mutsig_probabilities(df, mutsig_file):
    mutsig_df = pd.read_csv(mutsig_file)
    mutsig_df['mutsig_key'] = mutsig_df['Sample Names'] + '_' + mutsig_df['MutationTypes']
    mutsig_dict = dict(zip(mutsig_df['mutsig_key'], mutsig_df['mutsig']))
    df['mutsig'] = df['mutsig_key'].map(mutsig_dict)
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


def apply_step_5(file_path, assembly, mutsig_probabilities):
    
    df = pd.read_csv(file_path)
    
    # dtype optimization
    pred_columns = ['wt_prediction', 'mut_prediction']
    for col in pred_columns:
        df[col] = df[col].astype(np.float32)
    
    df = filter_rows_with_same_prediction(df)
    
    df = split_id_column(df)
    
    cat_columns = ['vcf_id','chr', 'pos', 'ref', 'alt', 'mirna_accession']
    for col in cat_columns:
        df[col] = df[col].astype('category')
    
    # column needed in pyensembl operations
    df["locus"] = df["chr"].astype(str) + ":" + df["pos"].astype(str)
    df = generate_gene_id_column(df, assembly)
    df['gene_id'] = df['gene_id'].astype('category')
    df = generate_is_intron_column(df, assembly)
    df.drop(columns=["locus"], inplace=True)
    
    # remove last character of vcf_id col if it is "+" to get the mutsigs
    df = generate_mutation_context_column(df)
    df.drop(columns=["chr", "pos"], inplace=True)
    
    df["mutation_context"] = df["mutation_context"].astype('category')
    
    
    df = add_mutsig_probabilities(df, mutsig_probabilities)
    df["mutsig"] = df["mutsig"].astype('category')
    df.drop(columns=["mutsig_key"], inplace=True)
    
    # deprecated for now
    # df["is_gain"] = df.mut_prediction > df.wt_prediction
    df["is_gene_upregulated"] = (df.wt_prediction > df.mut_prediction).astype('bool')

    # fill in columns that will produce nan
    for col in ["gene_id", "mutsig"]:
        existing_categories = df[col].cat.categories
        new_categories = existing_categories.append(pd.Index(['not_found']))
        df[col] = df[col].cat.set_categories(new_categories)
        df[col] = df[col].fillna('not_found')

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



def crawl_and_import_results_into_df(folder_path, ending_string, assembly):
    csv_files = []

    # Find CSV files
    for root, _, files in os.walk(folder_path):
        csv_files.extend(
            os.path.join(root, file)
            for file in files
            if file.endswith(f"{ending_string}.csv")
        )

    # Create a list to store DataFrames
    dataframes = []

    # Process each CSV file
    for file_path in csv_files:
        df = apply_step_5(file_path, assembly, MUTSIG_PROBABILITIES)
        dataframes.append(df)

    return pd.concat(dataframes, ignore_index=True)
