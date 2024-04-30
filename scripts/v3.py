import pandas as pd
from scripts import * 



def validate_ref_nucleotides_sharded(df, report_path, verbose=False):
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
        output_path (str, optional): The file path to write invalid rows to. Default is "invalid_rows.csv".
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


    if not invalid_rows.empty:
        if verbose:
            logging.warning(f"Writing {len(invalid_rows)} invalid rows to {report_path}")

        # Check if the file exists
        file_exists = os.path.isfile(report_path)

        # Open the file in append mode ('a') or write mode ('w')
        mode = 'a' if file_exists else 'w'
        with open(report_path, mode) as f:
            # If the file is new, write the header
            if not file_exists:
                f.write("id\n")

            # Write each invalid row to the file
            for _, row in invalid_rows.iterrows():
                f.write(f"{row['id']}\n")


    return df[~mask].drop("nuc_at_pos", axis=1)

def prepare_job_fastas_sharded(case_1, fasta_output_file):
    
    wild_type=True
    mirna_dict = pd.read_csv(MIRNA_CSV).set_index('mirna_accession')['sequence'].to_dict()
    
    sequence_column = 'wt_seq' if wild_type else 'mut_seq'
    mrna_dict = case_1.set_index('id')[sequence_column].to_dict()
    
    
    with open(fasta_output_file, 'w') as file:
        for string in generate_strings(mrna_dict, mirna_dict, wild_type):
            file.write(string)
        for string in generate_strings(mrna_dict, mirna_dict, not wild_type):
            file.write(string)


def run_rnaduplex_and_awk_sharded(input_file, output_file):

    try:
        with open(input_file, 'r') as input_f, tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
            subprocess.run(RNADUPLEX_LOCATION, stdin=input_f, stdout=temp_file, check=True, text=True)
            temp_file_path = temp_file.name
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.error(f"Error running RNAduplex on {input_file}: {str(e)}")
        return

    logging.debug(f"Running awk script on {input_file} results")

    try:
        awk_cmd = ["awk", "-f", AWK_SCRIPT_PATH, temp_file_path]
        awk_output = subprocess.check_output(awk_cmd, universal_newlines=True)
        with open(output_file, 'w') as output_f:
            output_f.write(awk_output)
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.error(f"Error running awk script on {temp_file_path}: {str(e)}")
    finally:
        try:
            os.remove(temp_file_path)
        except OSError as e:
            logging.warning(f"Error removing temporary file {temp_file_path}: {str(e)}")


def create_results_df(id_array, predictions, filter_range=QUANTILE_RANGE):
    df = pd.DataFrame({'id': id_array, 'prediction': predictions})
    df[['id', 'is_mutated']] = df['id'].str.rsplit('_', n=1, expand=True)
    df["is_mutated"] = df["is_mutated"].eq("mut")

    pivot_df = df.pivot(index='id', columns='is_mutated', values='prediction').reset_index()
    pivot_df.columns = ['id', 'wt_prediction', 'mut_prediction']
    pivot_df["pred_difference"] = (pivot_df["mut_prediction"] - pivot_df["wt_prediction"]).round(3)


    # filter out predictions that are outside of the quantile range
    mask = (pivot_df["pred_difference"].abs() >= filter_range)
    return pivot_df[mask]