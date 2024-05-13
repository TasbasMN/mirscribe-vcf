import os
from scripts.main_operations import *
from scripts.config import *
from scripts.cost_operations import *

def main():

    # Create the output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    processed_vcf = generate_csv_and_json_from_vcf(VCF_FULL_PATH, OUTPUT_DIR)
    
    print(processed_vcf)
    
    run_pipeline(processed_vcf, CHUNKSIZE, OUTPUT_DIR, VCF_ID)
    print("Processing completed. Results saved in the output directory.")

    results_filename = f"results_{VCF_ID}.csv"
    stitch_csv_files(OUTPUT_DIR, results_filename)
    remove_small_stitched_files(OUTPUT_DIR, results_filename)
    delete_fasta_files(OUTPUT_DIR)

    zip_file_path, zipped_files = zip_files(OUTPUT_DIR, "rnad", ".csv")
    delete_zipped_files(zip_file_path, zipped_files)


if __name__ == '__main__':
    main()
