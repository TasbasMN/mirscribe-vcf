import pandas as pd
from scripts.globals import MIRNA_CSV, RNADUPLEX_LOCATION
import subprocess, logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from multiprocessing import cpu_count


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')



def prepare_jobs_from_df(df):
    # Load the MIRNA_CSV file
    mirna_df = pd.read_csv(MIRNA_CSV)

    # Extract relevant data from the input DataFrame
    wt_sequences = df["wt_seq"].tolist()
    mutated_sequences = df["mut_seq"].tolist()
    identifiers = df["id"].tolist()
    mirna_identifiers = mirna_df["mirna_accession"].tolist()
    mirna_sequences = mirna_df["sequence"].tolist()

    # Create a dictionary to store miRNA sequences and identifiers
    mirna_dict = dict(zip(mirna_identifiers, mirna_sequences))

    # Create a generator for wild-type jobs
    def wt_job_generator():
        for wt, identifier in zip(wt_sequences, identifiers):
            for mirna_id, mirna in mirna_dict.items():
                yield (wt, mirna, identifier, mirna_id)

    # Create a generator for mutated jobs
    def mutated_job_generator():
        for mutated, identifier in zip(mutated_sequences, identifiers):
            for mirna_id, mirna in mirna_dict.items():
                yield (mutated, mirna, identifier, mirna_id)

    return wt_job_generator(), mutated_job_generator()


def run_rnaduplex_multithreaded(long_sequence, short_sequence, long_identifier, short_identifier):

    # handle the sequence input
    input_sequence = f"{long_sequence}\n{short_sequence}"

    result = subprocess.run(
        [RNADUPLEX_LOCATION, "-e", "5.0", "-s"],
        input=input_sequence,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,  # Capture stderr
        text=True)

    try:
        lines = result.stdout.split("\n")
        first_line = lines[0].split()
        brackets, long_start_end, _, short_start_end, energy = first_line
        long_bracket, short_bracket = brackets.split("&")
        start_long, end_long = map(int, long_start_end.split(","))
        start_short, end_short = map(int, short_start_end.split(","))
        energy_value = float(energy.strip("()"))

        return start_long - 1, end_long, long_bracket, start_short - 1, end_short, short_bracket, energy_value, long_identifier, short_identifier, long_sequence, short_sequence

    # predictions with positive energy prints out "( 5.0163)" with an extra whitespace in the beginning, so strip function adds another element "(" to the array.
    # we try to capture 6 elements with 5 assignments so the interpreter returns ValueError: too many values to unpack (expected 5). This part handles this issue.
    # predictions with negative energy prints out "(-5.0163)" so it works without problem.
    except ValueError as e:
        lines = result.stdout.split("\n")
        first_line = lines[0].split()

        if first_line == []:
            return 0, 1, ".", 0, 1, ".", 0, long_identifier, short_identifier, long_sequence, short_sequence
        
        brackets, long_start_end, _, short_start_end, _, energy = first_line
        long_bracket, short_bracket = brackets.split("&")
        start_long, end_long = map(int, long_start_end.split(","))
        start_short, end_short = map(int, short_start_end.split(","))
        energy_value = float(energy.strip("()"))

        return start_long - 1, end_long, long_bracket, start_short - 1, end_short, short_bracket, energy_value, long_identifier, short_identifier, long_sequence, short_sequence




def run_jobs_multithreaded(job_generator, binary_value, verbose=False):
    """
    Run jobs in a multithreaded manner using ProcessPoolExecutor.

    Args:
        job_generator (generator): A generator object that yields job tuples.
        binary_value (any): A binary value to be appended to the result of each job. 0 for wild-type, 1 for mutated.
        verbose (bool): If True, log the start and completion of each job.

    Returns:
        list: A list of tuples, where each tuple contains the result of a job and the binary value.
    """   
    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        futures = []
        for job in job_generator:
            future = executor.submit(run_rnaduplex_multithreaded, *job)
            futures.append((future, job))

        results = []
        for future, job in futures:
            if verbose:
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                logging.info(f"{timestamp} - Processing job: {job}")
            result = future.result()
            result_with_binary = result + (binary_value,)
            results.append(result_with_binary)
            if verbose:
                timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                logging.info(f"{timestamp} - Completed job: {job}")

    return results


def run_jobs_async(job_generator, binary_value, verbose=False, num_workers=cpu_count()):

    # Configure logging
    logging.basicConfig(level=logging.INFO if verbose else logging.WARNING)
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(run_rnaduplex_multithreaded, *job): job for job in job_generator}
        results = []

        for future in as_completed(futures):
            job = futures[future]
            try:
                result = future.result()
                result_with_binary = result + (binary_value,)
                results.append(result_with_binary)
                if verbose:
                    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                    logging.info(f"{timestamp} - Completed job: {job}")
                    
            except Exception as e:
                logging.error(f"Job {job} failed with error: {e}")
                continue

    return results