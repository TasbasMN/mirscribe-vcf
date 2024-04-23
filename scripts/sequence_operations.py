from scripts.globals import *


def calculate_au_content(sequence):
    au_count = sequence.count(
        'A') + sequence.count('T') + sequence.count('U')
    return None if len(sequence) == 0 else au_count / len(sequence)

def generate_mre_au_content_column(df):
    df["mre_au_content"] = df['mre_region'].apply(calculate_au_content)

    return df

def get_nucleotides_in_interval(chrom, start, end):
    """
    Given a chromosome name, start and end positions, this function reads the DNA sequence from the corresponding FASTA file and returns the nucleotides in the specified interval.

    Parameters:
    - chrom (str): The name of the chromosome.
    - start (int): The starting position of the interval.
    - end (int): The ending position of the interval.

    Returns:
    - nucleotides (str): The nucleotides in the specified interval.
    """
    file_path = f"{GRCH37_DIR}/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        start_offset = start - 1
        end_offset = end - 1
        num_start_new_lines = start_offset // line_length
        num_end_new_lines = end_offset // line_length
        start_byte_position = byte_position + start_offset + num_start_new_lines
        end_byte_position = byte_position + end_offset + num_end_new_lines
        file.seek(start_byte_position)

        # Read the nucleotides in the interval
        nucleotides = file.read(end_byte_position - start_byte_position + 1)

    # Remove newlines from the nucleotides
    nucleotides = nucleotides.replace('\n', '')

    return nucleotides


def get_nucleotide_at_position(chrom, position):
    """
    Given a chromosome name and a position, this function reads the DNA sequence from the corresponding FASTA file and returns the nucleotide at the specified position.

    Parameters:
    - chrom (str): The name of the chromosome.
    - position (int): The position of the nucleotide.

    Returns:
    - nucleotide (str): The nucleotide at the specified position.
    """
    file_path = f"{GRCH37_DIR}/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        offset = position - 1
        num_new_lines = offset // line_length
        byte_position = byte_position + offset + num_new_lines
        file.seek(byte_position)

        # Read the nucleotide at the position
        nucleotide = file.read(1)
    return nucleotide


def get_upstream_sequence(row, n=30):
    """
    Get the upstream sequence of length n from the given position.

    Args:
        row (pandas.Series): A row from the DataFrame containing the 'chr', 'pos', and 'ref_len' columns.
        n (int, optional): The length of the upstream sequence. Defaults to 30.

    Returns:
        str: The upstream sequence.
    """
    chrom = row['chr']
    pos = row['pos']
    upstream_start = max(1, pos - n)
    upstream_end = pos - 1
    return get_nucleotides_in_interval(chrom, upstream_start, upstream_end)


def get_downstream_sequence(row, n=30):
    """
    Get the downstream sequence of length n from the given position.

    Args:
        row (pandas.Series): A row from the DataFrame containing the 'chr', 'pos', and 'ref_len' columns.
        n (int, optional): The length of the downstream sequence. Defaults to 30.

    Returns:
        str: The downstream sequence.
    """
    chrom = row['chr']
    pos = row['pos']
    ref_len = len(row['ref'])
    downstream_start = pos + ref_len
    downstream_end = downstream_start + n - 1
    return get_nucleotides_in_interval(chrom, downstream_start, downstream_end)


def generate_mre_sequence(df):
    """
    Generate the miRNA response element (MRE) sequence for each row in the input DataFrame.

    Args:
        df (pandas.DataFrame): A DataFrame containing columns 'mrna_sequence', 'mrna_end', and 'mirna_start'.

    Returns:
        pandas.DataFrame: The input DataFrame with a new column 'mre_region' added.
    """
    # Calculate miRNA length
    df["mirna_length"] = df["mirna_sequence"].str.len()

    # Calculate MRE coordinates
    df["mre_end"] = df["mrna_end"] + df["mirna_start"]
    df["mre_start"] = df["mre_end"] - df["mirna_length"]

    # Ensure MRE start is not negative
    df["mre_start"] = df["mre_start"].apply(lambda x: max(x, 0))

    # Extract MRE sequence
    df["mre_region"] = df.apply(lambda row: row["mrna_sequence"][row["mre_start"]:row["mre_end"]], axis=1)

    # Drop temporary column
    df = df.drop(columns=["mirna_length"])

    return df


def add_sequence_columns(df):
    """
    Add columns for upstream, downstream, wild-type, and mutant sequences to a DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame containing 'chr', 'pos', 'ref', and 'alt' columns.

    Returns:
        pandas.DataFrame: The input DataFrame with additional columns for 'upstream_seq', 'downstream_seq',
                          'wt_seq', and 'mut_seq'.
    """

    df['upstream_seq'] = df.apply(lambda row: get_upstream_sequence(row, NUCLEOTIDE_OFFSET), axis=1)
    df['downstream_seq'] = df.apply(lambda row: get_downstream_sequence(row, NUCLEOTIDE_OFFSET), axis=1)
    df['wt_seq'] = df["upstream_seq"] + df["ref"] + df["downstream_seq"]
    df['mut_seq'] = df["upstream_seq"] + df["alt"] + df["downstream_seq"]

    return df