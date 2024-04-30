import argparse
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import os
from scripts import *


def process_chunk(chunk: pd.DataFrame, start_index: int, end_index: int, output_dir: str, vcf_id: str) -> tuple:
    start_time = datetime.now()
    # Perform your analysis pipeline on the chunk
    result = analysis_pipeline(chunk, start_index, end_index, output_dir, vcf_id)
    end_time = datetime.now()
    processing_time = end_time - start_time

    # Write the result to a CSV file in the output directory
    result_file = os.path.join(output_dir, f'result_{start_index}_{end_index}.csv')
    result.to_csv(result_file, index=False)

    return start_index, end_index, processing_time














def analysis_pipeline(df: pd.DataFrame, start_index: int, end_index: int, output_dir: str, vcf_id: str, verbose: bool = False,) -> pd.DataFrame:

    invalid_rows_report_file = os.path.join(output_dir, f"{vcf_id}_invalid_rows.csv")
    fasta_output_file = os.path.join(output_dir, f"{vcf_id}_{start_index}_{end_index}.fa")
    rnaduplex_output_file = os.path.join(output_dir, f"{vcf_id}_{start_index}_{end_index}.csv")

    df['id'] += '_' + df['chr'].astype(str) + '_' + df['pos'].astype(str) + '_' + df['ref'] + '_' + df['alt']
    df = validate_ref_nucleotides_sharded(df, invalid_rows_report_file)
    df = generate_is_mirna_column(df, grch=37)
    df = add_sequence_columns(df)
    case_1 = classify_and_get_case_1_mutations(df, vcf_id, start_index, end_index, output_dir)  
    prepare_job_fastas_sharded(case_1, fasta_output_file)
    run_rnaduplex_and_awk_sharded(fasta_output_file, rnaduplex_output_file)


    colnames = ["mutation_id", "mirna_accession", "mrna_dot_bracket_5to3", "mirna_dot_bracket_5to3", "mrna_start", "mrna_end", "mirna_start", "mirna_end", "pred_energy", "is_mutated"]
    df = pd.read_csv(rnaduplex_output_file, header=None, names=colnames)
    df["id"] = df["mutation_id"] + "_" + df["mirna_accession"] + "_" + df["is_mutated"]
    df = df.sort_values(by=['id', 'is_mutated'], ascending=[False, True])
    
    df = process_df_for_prediction(df)
    df, id_array = reorder_columns_for_prediction(df)
    predictions = make_predictions_with_xgb(df)
    df = create_results_df(id_array, predictions, filter_range=QUANTILE_RANGE)
    
    return df 




















def run_vcf_analysis(vcf_full_path: str, chunksize: int, output_dir: str, vcf_id: str):
    
    total_start_time = datetime.now()
    report_file = os.path.join(output_dir, 'analysis.csv')

    with ThreadPoolExecutor() as executor, open(report_file, 'w') as report:
        report.write(f"start_index,end_index,time_elapsed\n")
        report.write(f"start,time,{total_start_time}\n")
        

        futures = []
        start_index = 0
        for chunk in pd.read_csv(vcf_full_path, chunksize=chunksize, sep="\t", header=None, names=["chr", "pos", "id", "ref", "alt"]):
            end_index = start_index + len(chunk) - 1
            future = executor.submit(process_chunk, chunk, start_index, end_index, output_dir, vcf_id)
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
    parser = argparse.ArgumentParser(description='Process a VCF file in chunks using concurrent futures.')
    parser.add_argument('file_path', type=str, help='Path to the VCF file')
    parser.add_argument("-c", '--chunksize', default=100, type=int, help='Number of lines to process per chunk')
    parser.add_argument("-o", '--output_dir', type=str, default='./results', help='Path to the output directory')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')

    args = parser.parse_args()

    VCF_FULL_PATH = args.file_path
    VCF_ID = os.path.basename(VCF_FULL_PATH).split(".")[0]
    VERBOSE = args.verbose
    chunksize = args.chunksize
    
    OUTPUT_DIR = os.path.join(args.output_dir, f"{VCF_ID}_{chunksize}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)  # Create the output directory if it doesn't exist

    run_vcf_analysis(VCF_FULL_PATH, chunksize, OUTPUT_DIR, VCF_ID)
    print("Processing completed. Results saved in the output directory.")

if __name__ == '__main__':
    main()
