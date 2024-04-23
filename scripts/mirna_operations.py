import pandas as pd
from scripts.globals import MIRNA_COORDS_DIR, MIRNA_CSV
import os

def generate_is_mirna_column(df, grch):
    """
    Adds two columns to the input DataFrame:
    'is_mirna': 1 if the mutation falls within a miRNA region, 0 otherwise
    'mirna_accession': the accession number of the miRNA if 'is_mirna' is 1, None otherwise

    Args:
        df (pandas.DataFrame): The input DataFrame containing mutation data
        grch (int): The genome reference coordinate system version (e.g., 37, 38)

    Returns:
        pandas.DataFrame: The input DataFrame with two additional columns ('is_mirna' and 'mirna_accession')
    """
    # Construct the miRNA coordinates file path
    mirna_coords_file = os.path.join(MIRNA_COORDS_DIR, f"grch{grch}_coordinates.csv")

    # Load miRNA coordinates
    coords = pd.read_csv(mirna_coords_file)

    # Initialize new columns
    df['is_mirna'] = 0
    df['mirna_accession'] = None

    # Iterate over each mutation in the mutations dataframe
    for index, row in df.iterrows():
        mutation_chr = row['chr']
        mutation_start = row['pos']

        # Find matching miRNAs
        matching_rnas = coords.loc[(coords['chr'] == mutation_chr) &
                                   (coords['start'] <= mutation_start) &
                                   (coords['end'] >= mutation_start)]

        if not matching_rnas.empty:
            # Update the 'is_mirna' and 'mirna_accession' columns
            df.at[index, 'is_mirna'] = 1
            df.at[index, 'mirna_accession'] = matching_rnas['mirna_accession'].values[0]

    return df


def generate_mirna_conservation_column(df):
    """
    Add a 'mirna_conservation' column to the input DataFrame based on miRNA conservation data.

    Args:
        df (pandas.DataFrame): A DataFrame containing the 'mirna_accession' column.

    Returns:
        pandas.DataFrame: The input DataFrame with a 'mirna_conservation' column added.
    """
    mirna_df = (pd.read_csv(MIRNA_CSV, 
                            usecols=["mirna_accession", "conservation"])
                            .rename(columns={"conservation": "mirna_conservation"})
                            [["mirna_accession", "mirna_conservation"]])

    df = df.merge(mirna_df, on="mirna_accession", how="left")
    return df