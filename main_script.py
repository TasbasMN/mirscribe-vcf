import gc, os, logging
import pandas as pd

from scripts.dataframe_operations import *
from scripts.globals import *
from scripts.mirna_operations import *
from scripts.pyensembl_operations import *
from scripts.rnaduplex_operations import *
from scripts.sequence_operations import *
from scripts.truba_operations import *
from scripts.vcf_processing import *
from scripts.xgb_operations import *





def main():
    args = parse_cli_arguments()
    vcf_file_path = args.vcf
    start = args.start
    end = args.end
    output_dir = args.output_dir
    verbose = args.verbose
    
    # Extract the file name from the VCF file path
    vcf_id = os.path.basename(vcf_file_path).split(".")[0]

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    

    if verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    logging.info(f"VCF file: {vcf_file_path}")
    logging.info(f"Start index: {start}")
    logging.info(f"End index: {end}")
    logging.info(f"Output directory: {output_dir}")
    
    
    ## importing vcf and applying indices (if given)
    df = load_vcf_into_df(vcf_file_path, start, end)
    logging.info("VCF loaded into DataFrame.")
    print_memory_usage()

    
    ## modifying vcf dataframe
    augment_id(df)
    df = validate_ref_nucleotides(df)
    df = generate_is_mirna_column(df, grch=37)
    df = add_sequence_columns(df)
    generate_transcript_id_and_gene_name_columns(df, import_pyensembl(37))
    
    # classifying case 1 and case 2 mutations
    case_1 = df[df.is_mirna == 0][["id", "wt_seq", "mut_seq"]]
    case_2 = df[df.is_mirna == 1][["id", "wt_seq", "mut_seq"]]
    
    # saving case 1 and case 2 mutations to disk
    case_1.to_csv(os.path.join(output_dir, f"{vcf_id}_{start}_{end}_case_1.csv"), index=False)
    logging.info(f"Saved {len(case_1)} case 1 mutations to {os.path.join(output_dir, f'{vcf_id}_{start}_{end}_case_1.csv')}")
    
    if not case_2.empty:
        case_2.to_csv(os.path.join(output_dir, f"{vcf_id}_{start}_{end}_case_2.csv"), index=False)
        logging.info(f"Saved {len(case_2)} case 2 mutations to {os.path.join(output_dir, f'{vcf_id}_{start}_{end}_case_2.csv')}")
    else:
        logging.info(f"No case 2 mutations were found for {vcf_id}")
    
    
    # creating jobs from case 1 mutations and running them
    wt_jobs, mutated_jobs = prepare_jobs_from_df(df)
    wt_results = run_jobs_multithreaded(wt_jobs, 0, verbose)
    logging.info("Wild type jobs completed.")
    mutated_results = run_jobs_multithreaded(mutated_jobs, 1, verbose)
    logging.info("Mutated jobs completed.")
    
    # saving rnaduplex results to disk
    rnaduplex_results_file = os.path.join(output_dir, f"{vcf_id}_{start}_{end}_rnaduplex.csv")
    save_results_to_disk(wt_results, rnaduplex_results_file)
    save_results_to_disk(mutated_results, rnaduplex_results_file)
    logging.info("Wild type and mutated results saved to disk.")
    
    ############################# PREDICTION ################################
    
    # merging results into a df
    df = create_results_dataframe(wt_results, mutated_results)

    # adding features
    df = generate_alignment_string_from_dot_bracket(df)
    df = generate_match_count_columns(df)
    df = generate_ta_sps_columns(df)
    df = generate_mre_sequence(df)
    df = generate_important_sites(df)
    df = generate_mirna_conservation_column(df)
    df = generate_seed_type_columns(df)
    df = generate_mre_au_content_column(df)
    df = generate_local_au_content_column(df)
    logging.info("Features added to DataFrame.")
    
    # making predictions
    pred_df = make_predictions(df)
    logging.info("Predictions made.")

    # saving meaningful results
    meaningful_results_file = os.path.join(output_dir, f"{vcf_id}_{start}_{end}_meaningful_results.csv")
    pred_df[pred_df.pred_difference_binary != 0].to_csv(meaningful_results_file, index=False)
    logging.info("Meaningful results saved to disk.")


if __name__ == "__main__":
    main()
