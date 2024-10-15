import os
import csv
from typing import List
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

from scripts.pipeline_orchestration import *
from scripts.utils.misc_utils import time_it

from scripts.config import WORKERS, PROFILER


@time_it(enabled=PROFILER)
def run_pipeline(vcf_full_path: str, chunksize: int, output_dir: str, vcf_id: str):
    with ThreadPoolExecutor(max_workers=WORKERS) as executor:

        futures = []
        start_index = 0
        for chunk in pd.read_csv(vcf_full_path, chunksize=chunksize, sep="\t", header=None, names=["chr", "pos", "id", "ref", "alt"]):
            end_index = start_index + len(chunk) - 1
            future = executor.submit(
                process_chunk, chunk, start_index, end_index, output_dir, vcf_id)
            futures.append(future)
            start_index = end_index + 1

        for future in as_completed(futures):
            start_index, end_index = future.result()




def delete_fasta_files(directory: str):

    # List all files in the directory that end with .fa
    fasta_files = [f for f in os.listdir(directory) if f.endswith('.fa')]

    # Loop through the list of FASTA files and remove each one
    for file in fasta_files:
        file_path = os.path.join(directory, file)
        os.remove(file_path)


def delete_files(folder_path, prefix, extension):
    """
    Delete files in the specified folder that match both the given prefix and extension.

    Args:
    folder_path (str): Path to the folder containing files to delete.
    prefix (str): The prefix that matching files should start with.
    extension (str): The extension that matching files should end with.

    Returns:
    list: A list of filenames that were successfully deleted.
    """
    # Validate the folder path
    if not os.path.isdir(folder_path):
        raise ValueError(f"The folder path does not exist: {folder_path}")

    deleted_files = []

    # Iterate through files in the folder
    for filename in os.listdir(folder_path):
        # Check if the file matches both prefix and extension
        if filename.startswith(prefix) and filename.endswith(extension):
            file_path = os.path.join(folder_path, filename)
            
            # Attempt to delete the file
            try:
                os.remove(file_path)
                deleted_files.append(filename)
            except OSError as e:
                print(f"Error deleting file {filename}: {e}")


def stitch_and_cleanup_csv_files(output_dir: str, final_output_filename: str) -> None:
    """
    Stitch multiple result CSV files into one and remove the original files.
    
    Args:
        output_dir (str): Directory containing the CSV files.
        final_output_filename (str): Name of the final stitched file.
    """
    try:
        # Get and filter CSV files
        csv_files = [f for f in os.listdir(output_dir) 
                     if f.endswith('.csv') and f.startswith('result_')]
        csv_files.sort()  # Ensure consistent ordering
        
        if not csv_files:
            print("No matching CSV files found.")
            return

        final_output_path = os.path.join(output_dir, final_output_filename)
        removed_files: List[str] = []

        with open(final_output_path, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            header_written = False

            for filename in csv_files:
                file_path = os.path.join(output_dir, filename)
                with open(file_path, 'r', newline='') as infile:
                    reader = csv.reader(infile)
                    if not header_written:
                        header = next(reader)
                        writer.writerow(header)
                        header_written = True
                    else:
                        next(reader)  # Skip header in subsequent files
                    writer.writerows(reader)

    except Exception as e:
        print(f"An error occurred: {str(e)}")