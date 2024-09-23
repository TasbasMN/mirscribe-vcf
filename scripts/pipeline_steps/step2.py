import logging
import os
import pandas as pd
from scripts.globals import *
import tempfile
import subprocess


def classify_and_get_case_1_mutations(df, vcf_id, start, end, output_dir):
    """
    Classifies mutations into case 1 and case 2, saves case 2 mutations to disk,
    and returns the case 1 mutations.

    Args:
        df (pandas.DataFrame): DataFrame containing mutation data.
        vcf_id (str): ID of the VCF file.
        start (int): Start position of the region.
        end (int): End position of the region.
        output_dir (str): Path to the output directory.

    Returns:
        pandas.DataFrame: DataFrame containing case 1 mutations.
    """
    # Classify case 1 and case 2 mutations
    case_1 = df[df.is_mirna == 0][["id", "wt_seq", "mut_seq"]]
    case_2 = df[df.is_mirna == 1][["id", "wt_seq", "mut_seq"]]
    
    


    if not case_2.empty:
        # Save case 2 mutations to disk if any exist
        case_2_file = os.path.join(output_dir, f"{vcf_id}_{start}_{end}_case_2.csv")
        case_2.to_csv(case_2_file, index=False)
        print(f"case 2 mutations: {len(case_2)}")
        

    return case_1


def generate_fasta_representation_string(mrna_dict, mirna_dict, wild_type):
    """
    Generates FASTA type string based on the given mRNA and miRNA dictionaries and a flag indicating whether the sequence is wild type or mutant.

    Args:
        mrna_dict (dict): A dictionary containing mRNA sequences with their corresponding IDs as keys.
        mirna_dict (dict): A dictionary containing miRNA sequences with their corresponding IDs as keys.
        wild_type (bool): A flag indicating whether the sequence is wild type or mutant.

    Yields:
        str: A string in the format "{mRNA_ID}-{miRNA_ID}-wt" or "{mRNA_ID}-{miRNA_ID}-mut" followed by the mRNA and miRNA sequences separated by newline characters.

    """
    for mrna in mrna_dict:
        for mirna in mirna_dict:
            yield f">{mrna}-{mirna}-wt\n{mrna_dict[mrna]}\n{mirna_dict[mirna]}\n\n" if wild_type else f">{mrna}-{mirna}-mut\n{mrna_dict[mrna]}\n{mirna_dict[mirna]}\n\n"


def prepare_job_fastas_sharded(case_1, fasta_output_file):

    wild_type = True
    mirna_dict = pd.read_csv(MIRNA_CSV).set_index(
        'mirna_accession')['sequence'].to_dict()

    sequence_column = 'wt_seq' if wild_type else 'mut_seq'
    mrna_dict = case_1.set_index('id')[sequence_column].to_dict()

    with open(fasta_output_file, 'w') as file:
        for string in generate_fasta_representation_string(mrna_dict, mirna_dict, wild_type):
            file.write(string)
        for string in generate_fasta_representation_string(mrna_dict, mirna_dict, not wild_type):
            file.write(string)


def run_rnaduplex_and_awk_sharded(input_file, output_file):

    try:
        with open(input_file, 'r') as input_f, tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
            subprocess.run(RNADUPLEX_LOCATION, stdin=input_f,
                           stdout=temp_file, check=True, text=True)
            temp_file_path = temp_file.name
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.error(f"Error running RNAduplex on {input_file}: {str(e)}")
        return

    logging.debug(f"Running awk script on {input_file} results")

    try:
        awk_cmd = ["awk", "-f", AWK_SCRIPT_PATH, temp_file_path]
        awk_output = subprocess.check_output(awk_cmd, universal_newlines=True)
        with open(output_file, 'w') as output_f:
            output_f.write(awk_output)
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.error(
            f"Error running awk script on {temp_file_path}: {str(e)}")
    finally:
        try:
            os.remove(temp_file_path)
        except OSError as e:
            logging.warning(
                f"Error removing temporary file {temp_file_path}: {str(e)}")
