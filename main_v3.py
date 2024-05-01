import zipfile
import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import os
from scripts import *
import json


def time_it(func):
    """
    Decorator function that measures the execution time, memory usage, and system load average of a given function.

    Parameters:
        func (callable): The function to be timed.

    Returns:
        callable: The wrapped function that measures the execution time, memory usage, and system load average.
    """
    def wrapper(*args, **kwargs):
        process = psutil.Process()

        if hasattr(os, 'getloadavg'):
            load_before = os.getloadavg()  # System load average before function execution
        else:
            load_before = (None, None, None)

        start_time = time.time()  # Start timing
        memory_before = process.memory_info().rss / (1024 ** 2)

        exception = None  # Initialize exception to None
        try:
            result = func(*args, **kwargs)  # Call the function
        except Exception as exception:
            result = None

        end_time = time.time()  # End timing
        memory_after = process.memory_info().rss / (1024 ** 2)
        # Calculate duration in minutes
        duration = (end_time - start_time) / 60
        cpu_core_used = process.cpu_num()

        if hasattr(os, 'getloadavg'):
            load_after = os.getloadavg()  # System load average after function execution
        else:
            load_after = (None, None, None)

        # Prepare a concise representation of arguments for JSON output
        args_repr = [f"{len(a)}" if isinstance(
            a, pd.DataFrame) else repr(a) for a in args]

        # Create a dictionary to store function name, duration, and arguments
        log_entry = {
            "function_name": func.__name__,
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(start_time)),
            "duration_minutes": duration,
            "arguments": args_repr,
            "memory_before": memory_before,
            "memory_after": memory_after,
            "cpu_core_used": cpu_core_used,
            "exception": str(exception) if exception else None,
            "load_before": load_before,
            "load_after": load_after
        }

        # Write to JSON file
        json_file = os.path.join(OUTPUT_DIR, 'function_timings.json')
        with open(json_file, 'a') as f:
            json.dump(log_entry, f)
            f.write('\n')  # Write a newline to separate entries

        return result
    return wrapper


def zip_csv_files_with_vcf_id(directory: str, vcf_id: str, archive_name: str):
    """
    Compress all CSV files in the specified directory that start with the given vcf_id into a ZIP file.

    Parameters:
    - directory: The directory containing the CSV files to compress.
    - vcf_id: The specific prefix that the CSV files should start with to be included in the ZIP.
    - archive_name: The name of the resulting ZIP archive file.

    Returns:
    None
    """
    # Create a full path for the archive
    archive_path = os.path.join(directory, archive_name)

    # Create a new ZIP file
    with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        # Iterate over all the files in the directory
        for filename in os.listdir(directory):
            # Check if the file is a CSV and starts with the vcf_id
            if filename.endswith('.csv') and filename.startswith(vcf_id):
                # Create the full file path
                file_path = os.path.join(directory, filename)
                # Add file to the ZIP file
                zipf.write(file_path, filename)

    print(
        f"CSV files starting with {vcf_id} have been compressed into {archive_path}")


def stitch_csv_files(output_dir: str, final_output_filename: str):
    """
    Stitch together all CSV files in the specified output directory into a single CSV file
    without loading them all into memory at once.

    Parameters:
    - output_dir: The directory where the CSV files are stored.
    - final_output_filename: The name of the final stitched CSV file.

    Returns:
    None
    """
    # List all CSV files in the output directory that match the pattern
    csv_files = [f for f in os.listdir(
        output_dir) if f.endswith('.csv') and 'result_' in f]

    # Sort the files to maintain the order, if necessary
    csv_files.sort()

    # Define the path for the final output file
    final_output_path = os.path.join(output_dir, final_output_filename)

    # Open the final output file in write mode
    with open(final_output_path, 'w') as final_file:
        # Initialize a variable to track whether the header has been written
        header_written = False

        # Loop through the sorted CSV files
        for file in csv_files:
            file_path = os.path.join(output_dir, file)

            # Open the current CSV file in read mode
            with open(file_path, 'r') as read_file:
                # Read the header from the first file and write it to the final file
                header = read_file.readline()
                if not header_written:
                    final_file.write(header)
                    header_written = True

                # Stream the rest of the file's contents to the final file
                for line in read_file:
                    final_file.write(line)

    print(f"All CSV files have been stitched into {final_output_path}")


def remove_small_stitched_files(output_dir: str, exclude_filename: str):
    """
    Remove smaller stitched CSV files from the output directory, excluding the final stitched file.

    Parameters:
    - output_dir: The directory where the CSV files are stored.
    - exclude_filename: The name of the final stitched CSV file to exclude from deletion.

    Returns:
    None
    """
    # List all files in the output directory that match the pattern and are not the final stitched file
    files_to_remove = [f for f in os.listdir(output_dir) if f.startswith(
        'result_') and f.endswith('.csv') and f != exclude_filename]

    # Loop through the files and remove them
    for file in files_to_remove:
        file_path = os.path.join(output_dir, file)
        os.remove(file_path)

    print("Small stitched files removal complete.")


def delete_fasta_files(directory: str):
    """
    Delete all FASTA (.fa) files in the specified directory.

    Parameters:
    - directory: The directory to search for FASTA files.

    Returns:
    None
    """
    # List all files in the directory that end with .fa
    fasta_files = [f for f in os.listdir(directory) if f.endswith('.fa')]

    # Loop through the list of FASTA files and remove each one
    for file in fasta_files:
        file_path = os.path.join(directory, file)
        os.remove(file_path)

    print("All FASTA files have been deleted.")

def zip_files(folder_path, prefix, extension):
    """
    Creates a ZIP archive containing files in the specified folder
    that match the given prefix and extension.

    Args:
        folder_path (str): The path to the folder containing the files.
        prefix (str): The prefix to match in the file names.
        extension (str): The file extension to match (e.g., ".txt", ".py").

    Returns:
        tuple: A tuple containing the full path to the created ZIP archive file
               and a list of filenames that were zipped.
    """
    if not os.path.isdir(folder_path):
        raise ValueError(f"The folder path provided does not exist: {folder_path}")

    zip_filename = f"{prefix}_awk_outputs.zip"
    zip_filepath = os.path.join(folder_path, zip_filename)
    zipped_files = []

    with zipfile.ZipFile(zip_filepath, "w", compression=zipfile.ZIP_DEFLATED) as zip_archive:
        for file in os.listdir(folder_path):
            if file.startswith(prefix) and file.endswith(extension):
                file_path = os.path.join(folder_path, file)
                if os.path.isfile(file_path):
                    zip_archive.write(file_path, arcname=file)
                    zipped_files.append(file)

    if not zipped_files:
        os.remove(zip_filepath)
        return None, []
    
    return os.path.abspath(zip_filepath), zipped_files

def delete_zipped_files(zip_file_path, zipped_files):
    """
    Deletes the original files that were zipped, but only if the ZIP file
    was successfully created and is not empty.

    Args:
        zip_file_path (str): The full path to the created ZIP archive file.
        zipped_files (list): A list of filenames that were zipped.
    """
    if not os.path.isfile(zip_file_path) or os.path.getsize(zip_file_path) == 0:
        print(f"ZIP file not created or empty: {zip_file_path}")
        return

    folder_path = os.path.dirname(zip_file_path)
    for file in zipped_files:
        file_path = os.path.join(folder_path, file)
        if os.path.isfile(file_path):
            os.remove(file_path)
    
    print(f"Original files deleted after successful ZIP creation: {zip_file_path}")

@time_it
def process_chunk(chunk: pd.DataFrame, start_index: int, end_index: int, output_dir: str, vcf_id: str) -> tuple:
    start_time = datetime.now()
    # Perform your analysis pipeline on the chunk
    result = analysis_pipeline(
        chunk, start_index, end_index, output_dir, vcf_id)
    end_time = datetime.now()
    processing_time = end_time - start_time

    # Write the result to a CSV file in the output directory
    result_file = os.path.join(
        output_dir, f'result_{start_index}_{end_index}.csv')
    result.to_csv(result_file, index=False)

    return start_index, end_index, processing_time


def analysis_pipeline(df: pd.DataFrame, start_index: int, end_index: int, output_dir: str, vcf_id: str, verbose: bool = False,) -> pd.DataFrame:

    rnaduplex_output_file = os.path.join(
        output_dir, f"rnad_{vcf_id}_{start_index}_{end_index}.csv")

    if not SKIP_RNADUPLEX:
        invalid_rows_report_file = os.path.join(
            output_dir, f"invalid_rows_{vcf_id}.csv")
        fasta_output_file = os.path.join(
            output_dir, f"fasta_{vcf_id}_{start_index}_{end_index}.fa")

        df['id'] += '_' + df['chr'].astype(str) + '_' + df['pos'].astype(
            str) + '_' + df['ref'] + '_' + df['alt']
        df = validate_ref_nucleotides_sharded(df, invalid_rows_report_file)
        df = generate_is_mirna_column(df, grch=37)
        df = add_sequence_columns(df)
        case_1 = classify_and_get_case_1_mutations(
            df, vcf_id, start_index, end_index, output_dir)
        prepare_job_fastas_sharded(case_1, fasta_output_file)
        run_rnaduplex_and_awk_sharded(fasta_output_file, rnaduplex_output_file)

    colnames = ["mutation_id", "mirna_accession", "mrna_dot_bracket_5to3", "mirna_dot_bracket_5to3",
                "mrna_start", "mrna_end", "mirna_start", "mirna_end", "pred_energy", "is_mutated"]
    df = pd.read_csv(rnaduplex_output_file, header=None, names=colnames)
    df["id"] = df["mutation_id"] + "_" + \
        df["mirna_accession"] + "_" + df["is_mutated"]
    df = df.sort_values(by=['id', 'is_mutated'], ascending=[False, True])

    df = process_df_for_prediction(df)
    df, id_array = reorder_columns_for_prediction(df)
    predictions = make_predictions_with_xgb(df)
    df = create_results_df(id_array, predictions, filter_range=QUANTILE_RANGE)

    return df


@time_it
def run_vcf_analysis(vcf_full_path: str, chunksize: int, output_dir: str, vcf_id: str):

    # todo remove this timing code
    total_start_time = datetime.now()
    report_file = os.path.join(output_dir, 'analysis.csv')

    with ThreadPoolExecutor(max_workers=WORKERS) as executor, open(report_file, 'w') as report:
        report.write(f"start_index,end_index,time_elapsed\n")
        report.write(f"start,time,{total_start_time}\n")

        futures = []
        start_index = 0
        for chunk in pd.read_csv(vcf_full_path, chunksize=chunksize, sep="\t", header=None, names=["chr", "pos", "id", "ref", "alt"]):
            end_index = start_index + len(chunk) - 1
            future = executor.submit(
                process_chunk, chunk, start_index, end_index, output_dir, vcf_id)
            futures.append(future)
            start_index = end_index + 1

        for future in as_completed(futures):
            start_index, end_index, processing_time = future.result()
            report.write(f"{start_index},{end_index},{processing_time}\n")

    total_end_time = datetime.now()
    total_processing_time = total_end_time - total_start_time

    with open(report_file, 'a') as report:
        report.write(f"end,time,{total_end_time}\n")
        report.write(f"total,time: {total_processing_time}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Process a VCF file in chunks using concurrent futures.')
    parser.add_argument('file_path', type=str, help='Path to the VCF file')
    parser.add_argument("-c", '--chunksize', default=100,
                        type=int, help='Number of lines to process per chunk')
    parser.add_argument("-o", '--output_dir', type=str,
                        default='./results', help='Path to the output directory')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose logging')
    parser.add_argument('-w', '--workers', default=os.cpu_count(),
                        type=int, help='Number of concurrent workers')
    parser.add_argument('--skip-rnaduplex', action='store_true',
                        help='Skip RNAduplex analysis')

    args = parser.parse_args()

    VCF_FULL_PATH = args.file_path
    VCF_ID = os.path.basename(VCF_FULL_PATH).split(".")[0]
    VERBOSE = args.verbose
    chunksize = args.chunksize

    global WORKERS
    WORKERS = args.workers

    global SKIP_RNADUPLEX
    SKIP_RNADUPLEX = args.skip_rnaduplex

    global OUTPUT_DIR
    OUTPUT_DIR = os.path.join(args.output_dir, f"{VCF_ID}_{chunksize}")
    # Create the output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    run_vcf_analysis(VCF_FULL_PATH, chunksize, OUTPUT_DIR, VCF_ID)
    print("Processing completed. Results saved in the output directory.")

    results_filename = f"results_{VCF_ID}.csv"
    stitch_csv_files(OUTPUT_DIR, results_filename)
    remove_small_stitched_files(OUTPUT_DIR, results_filename)
    delete_fasta_files(OUTPUT_DIR)
    
    zip_file_path, zipped_files = zip_files(OUTPUT_DIR, "rnad", ".csv")
    delete_zipped_files(zip_file_path, zipped_files)

if __name__ == '__main__':
    main()
