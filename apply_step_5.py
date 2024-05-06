from scripts.pipeline_steps.step5 import *
import argparse
from scripts.globals import WT_THRESHOLD, MUT_THRESHOLD, MUTSIG_PROBABILITIES
import os
from pyensembl import EnsemblRelease

def find_csv_files(folder_path, starting_string):
    csv_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.startswith(starting_string):
                relative_path = os.path.relpath(os.path.join(root, file), folder_path)
                csv_files.append(relative_path)
    return csv_files

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find CSV files with a given ending string.")
    parser.add_argument("folder_path", help="Path to the folder to search")
    parser.add_argument("starting_string", help="Ending string for the CSV files")
    args = parser.parse_args()

    csv_files = find_csv_files(args.folder_path, args.starting_string)

    assembly = EnsemblRelease(75)



    for file_path in csv_files:
        file_path = f"results/{file_path}"
        print(file_path)


        df = apply_step_5(file_path, MUT_THRESHOLD, WT_THRESHOLD, MUTSIG_PROBABILITIES, assembly)

        csv_dir = os.path.dirname(file_path)
        VCF_ID = os.path.basename(file_path).split(".")[0]

        step5_csv = os.path.join(csv_dir, f"{VCF_ID}_step5.csv")
        df.to_csv(step5_csv, index=False)
