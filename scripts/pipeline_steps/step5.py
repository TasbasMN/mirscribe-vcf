import pandas as pd
from scripts.utils.sequence_utils import get_nucleotide_at_position
from scripts.pyensembl_operations import *
from scripts.globals import *


def filter_rows_by_thresholds(df, mut_threshold, wt_threshold):
    return df[(df["mut_prediction"] > mut_threshold) != (df["wt_prediction"] > wt_threshold)]


def split_id_column(df):
    df[['vcf_id', 'chr', 'pos', 'ref', 'alt', 'mirna_accession']
       ] = df['id'].str.split('_', expand=True)
    df.pos = df.pos.astype(int)
    return df


def generate_mutation_context_column(row):
    try:
        ref = row['ref']
        alt = row['alt']
        before = row['before']
        after = row['after']

        if ref in ['C', 'T']:  # Pyrimidine
            return f"{before}[{ref}>{alt}]{after}"
        elif ref in ['A', 'G']:  # Purine
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            ref_complement = complement[ref]
            alt_complement = complement[alt]
            before_complement = ''.join(
                complement[base] for base in before[::-1])
            after_complement = ''.join(
                complement[base] for base in after[::-1])
            return f"{after_complement}[{ref_complement}>{alt_complement}]{before_complement}"
        else:
            raise ValueError(f"Invalid reference nucleotide: {ref}")
    except KeyError as e:
        print(f"Error: Missing column - {str(e)}")
        return ''
    except Exception as e:
        print(f"Error: {str(e)}")
        return ''


def add_mutation_context_columns(df):
    df["before"] = df.apply(
        lambda x: get_nucleotide_at_position(x['chr'], x["pos"]-1), axis=1)
    df["after"] = df.apply(lambda x: get_nucleotide_at_position(
        x['chr'], x["pos"]+1), axis=1)
    df['mutation_context'] = df.apply(generate_mutation_context_column, axis=1)
    df["mutsig_key"] = df["vcf_id"] + "_" + df["mutation_context"]
    df.drop(columns=["before", "after", "ref", "alt"], inplace=True)
    return df


def add_mutsig_probabilities(df, mutsig_file):
    mutsig_df = pd.read_csv(mutsig_file)
    mutsig_df["mutsig_key"] = mutsig_df["Sample Names"] + \
        "_" + mutsig_df["MutationTypes"]
    mutsig_dict = mutsig_df.set_index('mutsig_key')['mutsig'].to_dict()
    df["mutsig"] = df["mutsig_key"].map(mutsig_dict)
    return df


def add_ensembl_data(df, ensembl_release):
    df['locus'] = df['chr'].astype(str) + ':' + df['pos'].astype(str)
    df = get_gene_ids(df, ensembl_release)
    df = get_gene_names(df, ensembl_release)
    df = get_transcript_ids(df, ensembl_release, canonical_only=True)
    df = get_gene_biotypes(df, ensembl_release)
    df.drop(columns=["chr", "pos", "locus"], inplace=True)
    return df


def merge_with_mirna_data(df):
    mirna_df = pd.read_csv(MIRNA_CSV, usecols=[
                           "mirna_accession", "mirna_name", "mirna_family"])
    df = df.merge(mirna_df, on="mirna_accession", how="left")
    return df


def apply_step_5(res_file, mut_threshold, wt_threshold, mutsig_file, ensembl_release):
    df = pd.read_csv(res_file)
    df = filter_rows_by_thresholds(df, mut_threshold, wt_threshold)
    df = split_id_column(df)
    df = add_mutation_context_columns(df)
    df = add_mutsig_probabilities(df, mutsig_file)
    df = add_ensembl_data(df, ensembl_release)
    df = merge_with_mirna_data(df)
    return df


def generate_is_intron_column(df, assembly):

    @lru_cache(maxsize=None)
    def cached_pyensembl_call(locus):

        chrom, pos = locus.split(':')
        result = assembly.exons_at_locus(chrom, int(pos))

        return np.nan if len(result) == 0 else result[0].exon_id

    unique_loci = df['locus'].unique()
    map_dict = {locus: cached_pyensembl_call(locus) for locus in unique_loci}
    df['is_exon'] = df['locus'].map(map_dict)
    
    mask = (df.is_exon.isna()) & (~df.gene_id.isna())

    df["is_intron"] = mask

    df.drop(columns=["is_exon"], inplace=True)

    return df