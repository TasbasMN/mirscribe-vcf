import numpy as np
import pandas as pd
import logging
import tempfile
import subprocess
import os
import concurrent.futures
import glob
import csv
import time

from scripts import *


def setup_logging_and_log_steps(output_dir, verbose=False):
    # Set up logging
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Define step descriptions
    step_descriptions = {
        "Step 1": "Load VCF into DataFrame",
        "Step 2": "Augment ID, validate ref nucleotides, generate is_mirna column, add sequence columns, classify and get case 1 mutations",
        "Step 3": "Prepare job FASTAs",
        "Step 4": "Process FASTA files in parallel and stitch CSV files",
        "Step 5": "Import combined CSV",
        "Step 6": "Generate feature columns",
        "Step 7": "Make predictions",
        "Step 8": "Save meaningful results to CSV"
    }

    def log_step_completion(step_name, step_start_time):
        end_time = time.time()
        duration = end_time - step_start_time
        logging.info(f"{step_name} completed in {duration:.2f} seconds")
        log_step_duration_to_file(step_name, duration)

    def log_step_duration_to_file(step_name, duration):
        description = step_descriptions.get(step_name, "")

        log_file = os.path.join(output_dir, "step_durations.csv")
        file_exists = os.path.isfile(log_file)
        with open(log_file, 'a', newline='') as csvfile:
            fieldnames = ['Step Name', 'Duration (seconds)', 'Description']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            if not file_exists:
                writer.writeheader()
            writer.writerow({'Step Name': step_name, 'Duration (seconds)': duration, 'Description': description})

    return log_step_completion

def classify_and_get_case_1_mutations(df, vcf_id, start, end, output_dir):
    """
    Classifies mutations into case 1 and case 2, saves case 2 mutations to disk,
    and returns the case 1 mutations.

    Args:
        df (pandas.DataFrame): DataFrame containing mutation data.
        vcf_id (str): ID of the VCF file.
        start (int): Start position of the region.
        end (int): End position of the region.
        output_dir (str): Path to the output directory.

    Returns:
        pandas.DataFrame: DataFrame containing case 1 mutations.
    """
    # Classify case 1 and case 2 mutations
    case_1 = df[df.is_mirna == 0][["id", "wt_seq", "mut_seq"]]
    case_2 = df[df.is_mirna == 1][["id", "wt_seq", "mut_seq"]]

    # Save case 2 mutations to disk if any exist
    case_2_file = os.path.join(output_dir, f"{vcf_id}_{start}_{end}_case_2.csv")
    if not case_2.empty:
        case_2.to_csv(case_2_file, index=False)
        logging.info(f"Saved {len(case_2)} case 2 mutations to {case_2_file}")
    else:
        logging.info(f"No case 2 mutations were found for {vcf_id}")

    logging.info(f"Found {len(case_1)} case 1 mutations.")

    return case_1




def generate_strings(mrna_dict, mirna_dict, wild_type):
    """
    Generates strings based on the given mRNA and miRNA dictionaries and a flag indicating whether the sequence is wild type or mutant.

    Args:
        mrna_dict (dict): A dictionary containing mRNA sequences with their corresponding IDs as keys.
        mirna_dict (dict): A dictionary containing miRNA sequences with their corresponding IDs as keys.
        wild_type (bool): A flag indicating whether the sequence is wild type or mutant.

    Yields:
        str: A string in the format "{mRNA_ID}-{miRNA_ID}-wt" or "{mRNA_ID}-{miRNA_ID}-mut" followed by the mRNA and miRNA sequences separated by newline characters.

    """
    for mrna in mrna_dict:
        for mirna in mirna_dict:
            yield f">{mrna}-{mirna}-wt\n{mrna_dict[mrna]}\n{mirna_dict[mirna]}\n\n" if wild_type else f">{mrna}-{mirna}-mut\n{mrna_dict[mrna]}\n{mirna_dict[mirna]}\n\n"


def generate_df_chunks(df, n_chunks):
    """
    Generates chunks of a DataFrame based on the specified number of chunks.

    Parameters:
        df (pandas.DataFrame): The DataFrame to be split into chunks.
        n_chunks (int): The number of chunks to generate.

    Yields:
        pandas.DataFrame: A chunk of the DataFrame.

    Examples:
        >>> df = pd.DataFrame({'A': [1, 2, 3, 4, 5], 'B': [6, 7, 8, 9, 10]})
        >>> for chunk in generate_df_chunks(df, 2):
        ...     print(chunk)
           A  B
        0  1  6
        1  2  7
           A  B
        2  3  8
        3  4  9
           A  B
        4  5  10
    """
    # Calculate the number of rows per chunk
    chunk_size = int(np.ceil(df.shape[0] / n_chunks))

    # Create an array of indices to split
    indices = np.arange(0, df.shape[0], chunk_size)

    # Use np.array_split to handle uneven division and yield chunks
    yield from np.array_split(df, indices[1:])


def prepare_job_fastas(df, output_dir, vcf_id):
    """
    Prepares FASTA files for a given DataFrame and writes them to the specified output directory.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data.
        output_dir (str): The directory where the FASTA files will be written.
        vcf_id (str): The ID of the VCF file.

    Returns:
        None

    Raises:
        Exception: If there is an error writing any of the FASTA files.

    Notes:
        - The function prepares FASTA files for both wild-type and mutated sequences.
        - The FASTA files are split into chunks based on the number of CPUs available.
        - Each chunk is written to a separate FASTA file with a unique name.
        - The function logs the successful writing of each chunk.
        - If there is an error writing any chunk, an error message is logged.
    """
    wild_type=True
    
    mirna_dict = pd.read_csv(MIRNA_CSV).set_index('mirna_accession')['sequence'].to_dict()

    for i, chunk in enumerate(generate_df_chunks(df, cpu_count())):
        sequence_column = 'wt_seq' if wild_type else 'mut_seq'
        mrna_dict = chunk.set_index('id')[sequence_column].to_dict()

        # Log the output file name and the indices
        start_index = chunk.index.min()
        end_index = chunk.index.max()
        
        output_file = os.path.join(output_dir, f"{vcf_id}_chunk_{i}_{start_index}_{end_index}.fa")

        try:
            with open(output_file, 'w') as file:
                for string in generate_strings(mrna_dict, mirna_dict, wild_type):
                    file.write(string)
                for string in generate_strings(mrna_dict, mirna_dict, not wild_type):
                    file.write(string)
            logging.debug(f'Wrote chunk {i}, indices {start_index}-{end_index}')
        except Exception as e:
            logging.error(f'Error writing chunk {i}, indices {start_index}-{end_index}: {e}')
                
    
def run_rnaduplex_and_awk(input_file, awk_script_path, output_file):
    """
    Run RNAduplex on the input file and then run an awk script on the results.
    
    Args:
        input_file (str): The path to the input file.
        awk_script_path (str): The path to the awk script.
        output_file (str): The path to the output file.
    
    Returns:
        None: If there is an error running RNAduplex or awk script.
    
    Raises:
        FileNotFoundError: If RNAduplex or awk script is not found.
        subprocess.CalledProcessError: If there is an error running RNAduplex or awk script.
    
    """
    logging.debug(f"Running RNAduplex on {input_file}")

    try:
        with open(input_file, 'r') as input_f, tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
            subprocess.run(RNADUPLEX_LOCATION, stdin=input_f, stdout=temp_file, check=True, text=True)
            temp_file_path = temp_file.name
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.error(f"Error running RNAduplex on {input_file}: {str(e)}")
        return

    logging.debug(f"Running awk script on {input_file} results")

    try:
        awk_cmd = ["awk", "-f", awk_script_path, temp_file_path]
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


def process_fasta_files_parallel(work_dir):
    """
    Process fasta files in parallel using ProcessPoolExecutor.
    
    Parameters:
    work_dir (str): The directory containing the fasta files to be processed.
    
    Returns:
    None
    """
    input_files = [os.path.join(work_dir, file) for file in os.listdir(work_dir) if file.endswith(".fa")]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for input_file in input_files:
            output_file = os.path.join(work_dir, f"{os.path.basename(input_file)}.csv")
            future = executor.submit(run_rnaduplex_and_awk, input_file, AWK_SCRIPT_PATH, output_file)
            futures.append(future)

        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logging.error(f"Error processing file: {str(e)}")


def stitch_csv_files(work_dir, vcf_id, start, end):
    logging.debug(f"Stitching CSV files from {work_dir}")

    # Get a sorted list of CSV files in the output folder
    csv_files = sorted(glob.glob(os.path.join(work_dir, f"{vcf_id}*.csv")))

    if not csv_files:
        logging.warning(f"No CSV files starting with '{vcf_id}' found in {work_dir}")
        return

    logging.debug(f"Found {len(csv_files)} CSV files starting with '{vcf_id}'")

    # Check if the combined file exists, and if so, clear it
    combined_csv = os.path.join(work_dir, f"combined_{vcf_id}_{start}_{end}.csv")
    
    if os.path.exists(combined_csv):
        logging.info(f"Clearing existing file: {combined_csv}")
        with open(combined_csv, "w", newline="") as f:
            pass

    # Open the combined CSV file in write mode
    with open(combined_csv, "a", newline="") as f:
        writer = csv.writer(f)

        # Iterate over each CSV file
        for csv_file in csv_files:
            logging.debug(f"Writing {csv_file} to {combined_csv}")

            with open(csv_file, "r") as input_csv:
                reader = csv.reader(input_csv)

                # Write each row from the input CSV to the combined CSV
                for row in reader:
                    writer.writerow(row)

            # Delete the stitched CSV fragment
            os.remove(csv_file)
            logging.debug(f"Deleted {csv_file}")

    logging.debug(f"Combined CSV file created: {combined_csv}")
    
    return combined_csv


def import_combined_csv(combined_csv):
    logging.debug(f"Processing combined csv from {combined_csv}")

    colnames = ["mutation_id", "mirna_accession", "mrna_dot_bracket_5to3", "mirna_dot_bracket_5to3", "mrna_start", "mrna_end", "mirna_start", "mirna_end", "pred_energy", "is_mutated"]
    try:
        df = pd.read_csv(combined_csv, header=None, names=colnames)
    except FileNotFoundError:
        logging.error(f"File not found: {combined_csv}")
        return None

    # augments id
    df["id"] = df["mutation_id"] + "_" + df["mirna_accession"] + "_" + df["is_mutated"]
    df = df.sort_values("id")

    logging.debug("Combined csv imported")
    return df


def downcast_df(df):
    # Log the initial RAM allocation
    logging.info(f"Initial RAM allocation: {df.memory_usage(deep=True).sum() / (1024 ** 2):.2f} MB")

    for column in df.columns:
        if df[column].dtype == 'object':
            df[column] = df[column].astype('category')
        elif df[column].dtype == 'float64':
            df[column] = pd.to_numeric(df[column], downcast='float')
        else:
            df[column] = pd.to_numeric(df[column], downcast='integer')

    # Log the final RAM allocation
    logging.info(f"Final RAM allocation: {df.memory_usage(deep=True).sum() / (1024 ** 2):.2f} MB")

    return df




import pandas as pd

def process_df_for_prediction(df):
    """
    Process the input dataframe by adding various columns and performing operations.

    Args:
        df (pandas.DataFrame): The input dataframe.

    Returns:
        pandas.DataFrame: The processed dataframe.
    """
    ## Add miRNA sequence
    mirna_dict = pd.read_csv(MIRNA_CSV).set_index('mirna_accession')['sequence'].to_dict()
    df["mirna_sequence"] = df["mirna_accession"].map(mirna_dict)

    ## Add miRNA conservation
    df = generate_mirna_conservation_column(df)
    df.drop("mirna_accession", axis=1, inplace=True)

    ## Add mRNA sequence
    df = split_mutation_ids(df)
    df['is_mutated'] = df['is_mutated'].isin(['mt', 'mut'])
    df = add_sequence_columns(df)
    df['mrna_sequence'] = df['wt_seq'].where(~df['is_mutated'], df['mut_seq'])
    column_names = ["chr", "pos", "ref", "alt", "upstream_seq", "downstream_seq", "wt_seq", "mut_seq", "is_mutated"]
    df.drop(columns=column_names, inplace=True)

    ## Generate MRE sequence
    df = generate_mre_sequence(df)
    df.drop(columns=["mrna_start", "mrna_end", "mre_start", "mre_end"], inplace=True)

    ## Generate local AU content
    df = generate_local_au_content_column(df)
    df.drop("mrna_sequence", axis=1, inplace=True)

    ## Generate TA SPS columns
    df = generate_ta_sps_columns(df)

    ## Generate alignment string from dot bracket
    df = generate_alignment_string_from_dot_bracket(df)
    df.drop(columns=["mirna_start", "mirna_end", "mirna_sequence"], inplace=True)

    ## Generate match count columns and important sites
    df = generate_match_count_columns(df)
    df = generate_important_sites_optimized(df)
    df = generate_seed_type_columns(df)
    df.drop("alignment_string", axis=1, inplace=True)

    ## Generate MRE AU content
    df = generate_mre_au_content_column(df)
    df.drop("mre_region", axis=1, inplace=True)

    return df
