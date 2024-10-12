import psutil
import os
import zipfile
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


def zip_files(folder_path, prefix, extension):

    if not os.path.isdir(folder_path):
        raise ValueError(
            f"The folder path provided does not exist: {folder_path}")

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

    print(f"Files zipped successfully: {zip_filepath}")
    return os.path.abspath(zip_filepath), zipped_files


def delete_zipped_files(zip_file_path, zipped_files):

    if not os.path.isfile(zip_file_path) or os.path.getsize(zip_file_path) == 0:
        print(f"ZIP file not created or empty: {zip_file_path}")
        return

    folder_path = os.path.dirname(zip_file_path)
    for file in zipped_files:
        file_path = os.path.join(folder_path, file)
        if os.path.isfile(file_path):
            os.remove(file_path)

    print(
        f"Original files deleted after successful ZIP creation: {zip_file_path}")


def print_memory_usage(verbose=False):

    if verbose:
        # Get the current process
        process = psutil.Process()

        # Get the memory info for the current process
        mem_info = process.memory_info()

        # Get the total memory and used memory in bytes
        total_memory = psutil.virtual_memory().total
        used_memory = mem_info.rss

        # Calculate the memory usage percentage
        memory_percentage = (used_memory / total_memory) * 100

        # Log the memory usage information
        logging.info(
            f"Used Memory: {used_memory / 1024 / 1024:.2f} MB, memory percentage: {memory_percentage:.2f}%")


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