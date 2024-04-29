from scripts import *
import argparse
import cProfile
import shutil
import time
import csv
import os
import sys


def main():
    pipeline_start_time = time.time()
    
    
    # Step 1
    step_start_time = time.time()
    df = load_vcf_into_df(vcf_file_path, start, end)
    log_step_completion("Step 1", step_start_time)


    # Step 2
    step_start_time = time.time()
    augment_id(df)
    df = validate_ref_nucleotides(df)
    df = generate_is_mirna_column(df, grch=37)
    df = add_sequence_columns(df)
    case_1 = classify_and_get_case_1_mutations(df, vcf_id, start, end, output_dir)
    log_step_completion("Step 2", step_start_time)


    # Step 3
    step_start_time = time.time()
    prepare_job_fastas(case_1, output_dir, vcf_id)
    log_step_completion( "Step 3", step_start_time)
    
    
    # Step 4
    step_start_time = time.time()
    process_fasta_files_parallel(output_dir)
    combined_csv = stitch_csv_files(output_dir, vcf_id, start, end)
    log_step_completion("Step 4", step_start_time)


    # Step 5
    step_start_time = time.time()
    df = import_combined_csv(combined_csv)
    log_step_completion("Step 5", step_start_time)
    

    # Step 6
    step_start_time = time.time()
    df = process_df_for_prediction(df)
    log_step_completion("Step 6", step_start_time)
    
    
   # Step 7
    step_start_time = time.time()
    df, id_array = reorder_columns_for_prediction(df)
    predictions = make_predictions_with_xgb(df)
    pred_df = create_prediction_result_df(id_array, predictions, filter_range=QUANTILE_RANGE, brief_output=False)
    log_step_completion("Step 7", step_start_time)
    

   # Step 8
    step_start_time = time.time()
    meaningful_results_file = os.path.join(output_dir, f"{vcf_id}_{start}_{end}_meaningful_results.csv")
    pred_df.to_csv(meaningful_results_file, index=False)
    log_step_completion("Step 8", step_start_time)
 
    
    # Log the total duration of the job
    end_time = time.time()
    total_duration = end_time - pipeline_start_time
    logging.info(f"Total job duration: {total_duration:.2f} seconds")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="VCF file to analyze")
    parser.add_argument("-s", "--start", default=0, type=int, help="start index")
    parser.add_argument("-e", "--end", default=-1, type=int, help="end index")
    parser.add_argument("-o", "--output_dir", default="./results", help="Output directory for runtime file, defaults to ./results")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("-p", "--profile", action='store_true', help='Enable profiler')

    args = parser.parse_args()
    vcf_file_path = args.vcf
    start = args.start
    end = args.end
    profile = args.profile
    
    # Extract the file name from the VCF file path
    vcf_id = os.path.basename(vcf_file_path).split(".")[0]
    output_dir = os.path.join(args.output_dir, f"{vcf_id}_{start}_{end}")
    
    log_step_completion = setup_logging_and_log_steps(output_dir, args.verbose)
    

    # Remove the output directory if it already exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # Create a new output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # profiler configuration
    if args.profile:
        logging.info("Enabling profiler.")
        profiler = cProfile.Profile()
        profiler.enable()
    
    main()
    
    if args.profile:
        profiler.disable()
        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        
        profiler_dump = f"profile_output_{timestamp}.pstats"
        profiler.dump_stats(os.path.join(output_dir, profiler_dump))
        logging.info(f"Profiler results saved to {os.path.join(output_dir, profiler_dump)}")
        
        