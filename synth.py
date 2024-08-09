import os
from scripts.main_operations import *
from scripts.config import *
import time
from memory_profiler import profile

@profile
def main():

    # Create the output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # print and save current time to a variable
    start = time.time()
    
    
    run_pipeline(VCF_FULL_PATH, CHUNKSIZE, OUTPUT_DIR, VCF_ID)
    print("Processing completed. Results saved in the output directory.")

    results_filename = f"results_{VCF_ID}.csv"
    stitch_csv_files(OUTPUT_DIR, results_filename)
    remove_small_stitched_files(OUTPUT_DIR, results_filename)
    delete_fasta_files(OUTPUT_DIR)

    zip_file_path, zipped_files = zip_files(OUTPUT_DIR, "rnad", ".csv")
    delete_zipped_files(zip_file_path, zipped_files)

    end = time.time()
    print(f"Total runtime: {end - start:.2f} seconds")
    
if __name__ == '__main__':
    main()
