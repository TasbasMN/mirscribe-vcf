import pandas as pd
import numpy as np
from scripts.sequence_operations import get_nucleotide_at_position, get_nucleotides_in_interval

def load_vcf_into_df(file_path, start=0, end=-1):
    """
    Load a VCF (Variant Call Format) file into a pandas DataFrame.
    
    Parameters:
    - file_path (str): The full path to the VCF file to be loaded.
    - start (int): The starting index for slicing the DataFrame (default: 0).
    - end (int): The ending index for slicing the DataFrame (default: -1, which means no end limit).
    
    Returns:
    - DataFrame: A pandas DataFrame containing the VCF file data, with columns for 'chr', 'pos', 'id', 'ref', and 'alt', sliced based on the start and end indices.
    
    Raises:
    - FileNotFoundError: If the specified file does not exist.
    - ValueError: If the file is not in the expected VCF format.
    - Exception: For other unforeseen errors during file loading.
    """
    # Ensure the file path ends with '.vcf', indicating it is a VCF file
    if not file_path.endswith('.vcf'):
        raise ValueError("The file specified does not appear to be a VCF file. Please provide a valid .vcf file.")

    try:
        df = pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            names=["chr", "pos", "id", "ref", "alt"],
            comment='#',
        )

        if end == -1:
            end = len(df)

        return df[start:end]

    except FileNotFoundError as e:
        # Raise a more descriptive error if the file is not found
        raise FileNotFoundError(f"The file at {file_path} was not found.") from e
    
    except pd.errors.EmptyDataError as e:
        # Handle the case where the file is empty or does not contain expected data
        raise ValueError(
            f"The file at {file_path} is empty or not in the expected VCF format."
        ) from e
    except Exception as e:
        # Catch-all for other exceptions, with a generic error message
        print(f"An error occurred while trying to read the file at {file_path}: {e}")
        raise
    
    
def augment_id(df):
    """Augment the 'id' column with additional information."""
    df['id'] += '_' + df['chr'].astype(str) + '_' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt']


import logging

def validate_ref_nucleotides(df, output_file="invalid_rows.csv", verbose=False):
    """
    Check if the 'ref' column matches the 'nucleotide_at_position' column, fetched from FASTA, in a DataFrame.
    Write invalid rows to a file and return the valid rows for downstream analysis.

    The 'nucleotide_at_position' column is populated as follows:
    - For rows where the reference length (ref_len) is greater than 1, the get_nucleotides_in_interval function
      is used to fetch the nucleotides in the interval [pos, pos + ref_len - 1] for the given chromosome (chr).
    - For rows where the reference length (ref_len) is 1, the get_nucleotide_at_position function is used to
      fetch the nucleotide at the given position (pos) for the given chromosome (chr).
    - For all other cases, an empty string is assigned.

    Args:
        df (pandas.DataFrame): The input DataFrame.
        output_file (str, optional): The file name to write invalid rows to. Default is "invalid_rows.csv".
        verbose (bool, optional): If True, log messages indicating the progress. Default is False.

    Returns:
        pandas.DataFrame: The DataFrame containing valid rows with matching 'ref' and 'nucleotide_at_position',
                           without the 'nucleotide_at_position' column.
    """
    if verbose:
        logging.info("Validating reference nucleotides...")

    # Add ref_len and alt_len columns
    df["ref_len"] = df["ref"].str.len()
    df["alt_len"] = df["alt"].str.len()

    # For rows where the reference length (ref_len) is greater than 1:
    # - Use the get_nucleotides_in_interval function to fetch the nucleotides
    #   in the interval [pos, pos + ref_len - 1] for the given chromosome (chr)
    df['nuc_at_pos'] = np.where(
        df['ref_len'] > 1,
        df.apply(lambda x: get_nucleotides_in_interval(x['chr'], x['pos'], x["pos"] + x["ref_len"] - 1), axis=1),
        # For rows where the reference length (ref_len) is 1:
        # - Use the get_nucleotide_at_position function to fetch the nucleotide
        #   at the given position (pos) for the given chromosome (chr)
        np.where(
            df['ref_len'] == 1,
            df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']), axis=1),
            # For all other cases, set the value to an empty string
            ""
        )
    )

    ## Check if ref matches nucleotide_at_position
    mask = df['ref'] != df['nuc_at_pos']

    ## Isolate invalid rows
    invalid_rows = df[mask]

    ## Write invalid rows to a file
    if not invalid_rows.empty:
        if verbose:
            logging.warning(f"Writing {len(invalid_rows)} invalid rows to {output_file}")
        invalid_rows[["id"]].to_csv(output_file, index=False)

    if verbose:
        logging.info(f"Found {len(df[~mask])} valid rows.")

    return df[~mask].drop("nuc_at_pos", axis=1)
