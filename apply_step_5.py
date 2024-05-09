# Import necessary modules and functions
from scripts.pipeline_steps.step5 import *  # Import all functions from step5
import argparse
from scripts.globals import WT_THRESHOLD, MUT_THRESHOLD, MUTSIG_PROBABILITIES
import os
from pyensembl import EnsemblRelease

# Function to find CSV files with a given starting string
def find_csv_files(folder_path, starting_string):
    """
    Find CSV files with a given starting string in a folder.

    Args:
        folder_path (str): Path to the folder to search
        starting_string (str): Starting string for the CSV files

    Returns:
        list: List of absolute paths to the CSV files
    """
    csv_files = []
    for root, _, files in os.walk(folder_path):
        csv_files.extend(
            os.path.join(root, file)
            for file in files
            if file.startswith(starting_string)
        )
    return csv_files

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Find CSV files with a given starting string.")
    parser.add_argument("folder_path", help="Path to the folder to search")
    parser.add_argument("starting_string", help="Ending string for the CSV files")
    args = parser.parse_args()

    # Find CSV files
    csv_files = find_csv_files(args.folder_path, args.starting_string)

    # Initialize EnsemblRelease
    assembly = EnsemblRelease(75)

    # Process each CSV file
    for file_path in csv_files:
        print(file_path)

        # Apply step 5 to the CSV file
        df = apply_step_5(file_path, MUT_THRESHOLD, WT_THRESHOLD, MUTSIG_PROBABILITIES, assembly)

        # Save the result to a new CSV file
        csv_dir = os.path.dirname(file_path)
        VCF_ID = os.path.basename(file_path).split(".")[0]
        step5_csv = os.path.join(csv_dir, f"{VCF_ID}_step5.csv")
        df.to_csv(step5_csv, index=False)