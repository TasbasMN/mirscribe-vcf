import pandas as pd
import os
import gc
from scripts.utils.misc_utils import time_it

from scripts.pipeline_steps.step1 import *
from scripts.pipeline_steps.step2 import *
from scripts.pipeline_steps.step3 import *
from scripts.pipeline_steps.step4 import *

from scripts.config import SKIP_RNADUPLEX, FILTER_THRESHOLD, PROFILER

@time_it(enabled=PROFILER)
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

        # Step 1: Data Preprocessing
        df = validate_ref_nucleotides_sharded(df, invalid_rows_report_file)
        df = generate_is_mirna_column(df, grch=37)
        df = add_sequence_columns(df)

        # Step 2: Data Processing
        case_1 = classify_and_get_case_1_mutations(
            df, vcf_id, start_index, end_index, output_dir)
        
        # gc
        del df
        gc.collect()
        
        prepare_job_fastas_sharded(case_1, fasta_output_file)
        
        # gc
        del case_1
        gc.collect()
        
        run_rnaduplex_and_awk_sharded(fasta_output_file, rnaduplex_output_file)

        # Step 3: Prediction Preprocessing
        df = process_rnaduplex_output(rnaduplex_output_file)
        df = generate_mirna_conservation_column(df)
        df.drop("mirna_accession", axis=1, inplace=True)
        df = split_mutation_ids(df)
        df['is_mutated'] = df['is_mutated'].isin(['mt', 'mut'])
        df = add_sequence_columns(df)
        df['mrna_sequence'] = df['wt_seq'].where(
            ~df['is_mutated'], df['mut_seq'])
        column_names = ["chr", "pos", "ref", "alt", "upstream_seq",
                        "downstream_seq", "wt_seq", "mut_seq", "is_mutated"]
        df.drop(columns=column_names, inplace=True)
        
        
        df = generate_mre_sequence_column(df)
        
        # add a mask that checks if the mutation is in the MRE region
        mask_if_mutation_in_mre = (df.mrna_start < 32) & (df.mrna_end > 30)
        df["is_mutation_in_mre"] = mask_if_mutation_in_mre
        df.drop(columns=["mrna_start", "mrna_end",
                "mre_start", "mre_end"], inplace=True)
        
        df = generate_local_au_content_column(df)
        df.drop("mrna_sequence", axis=1, inplace=True)
        
        df = generate_ta_sps_columns(df)
        df = generate_alignment_string_from_dot_bracket(df)
        df.drop(columns=["mirna_start", "mirna_end",
                "mirna_sequence"], inplace=True)
        
        df = generate_match_count_columns(df)
        df = generate_important_sites_column(df)
        df = generate_seed_type_columns(df)
        df.drop("alignment_string", axis=1, inplace=True)
        
        df = generate_mre_au_content_column(df)
        df.drop("mre_region", axis=1, inplace=True)

        # Step 4: Prediction
        df, id_array, binary_array = reorder_columns_for_prediction(df)

        predictions = make_predictions_with_xgb(df)
        
        # gc
        del df
        gc.collect()
        
        df = create_results_df(id_array, predictions, binary_array,
                               filter_range=FILTER_THRESHOLD)
        
        df.drop(columns=["binary_array"], inplace=True)

    return df


# @time_it(enabled=PROFILER)
def process_chunk(chunk: pd.DataFrame, start_index: int, end_index: int, output_dir: str, vcf_id: str) -> tuple:

    result = analysis_pipeline(
        chunk, start_index, end_index, output_dir, vcf_id)

    # Write the result to a CSV file in the output directory
    result_file = os.path.join(
        output_dir, f'result_{start_index}_{end_index}.csv')
    result.to_csv(result_file, index=False)

    return start_index, end_index



