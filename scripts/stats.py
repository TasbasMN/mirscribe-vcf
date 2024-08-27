import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact


def safe_log2_ratio(df):
    return np.log2((df[True] + 1) / (df[False] + 1))


def process_data(df, name):
    grouped = df.groupby(['gene_name', 'is_gene_upregulated']
                         ).size().unstack(fill_value=0)

    grouped['log2_ratio'] = safe_log2_ratio(grouped)

    return grouped.sort_values('log2_ratio', ascending=False)


def compare_datasets(real_df, synth_df):
    common_genes = list(set(real_df.index) & set(synth_df.index))
    print(f"\nNumber of common genes: {len(common_genes)}")

    real_subset = real_df.loc[common_genes]
    synth_subset = synth_df.loc[common_genes]

    result = pd.DataFrame({
        'real_log2_ratio': real_subset['log2_ratio'],
        'synth_log2_ratio': synth_subset['log2_ratio'],
        'real_upregulated': real_subset[True],
        'real_downregulated': real_subset[False],
        'synth_upregulated': synth_subset[True],
        'synth_downregulated': synth_subset[False]
    })

    result['log2fc'] = result['real_log2_ratio'] - \
        result['synth_log2_ratio']

    # Calculate z-scores
    epsilon = 1e-10
    result['z_score'] = (result['log2fc'] - result['log2fc'].mean()
                         ) / (result['log2fc'].std() + epsilon)

    # Calculate Fisher's exact test p-values
    result['fisher_p_value'] = result.apply(lambda row: fisher_exact([[row['real_upregulated'], row['real_downregulated']],
                                                               [row['synth_upregulated'], row['synth_downregulated']]])[1], axis=1)

    # Adjust p-values for multiple testing
    result['adjusted_fisher_p_value'] = multipletests(
        result['fisher_p_value'], method='fdr_bh')[1]

    print(f"\nShape of result dataframe: {result.shape}")
    print(
        f"Range of log2fc: {result['log2fc'].min()} to {result['log2fc'].max()}")
    print(
        f"Range of z_scores: {result['z_score'].min()} to {result['z_score'].max()}")

    return result


def filter_and_categorize(df, z_threshold=2, log2_threshold=0.32, p_threshold=0.05):
    filtered = df[(abs(df['z_score']) > z_threshold) &
                  (abs(df['log2fc']) > log2_threshold) &
                  (df['adjusted_fisher_p_value'] < p_threshold)].copy()
    print(f"\nNumber of genes after filtering: {len(filtered)}")

    if len(filtered) == 0:
        print("No genes pass the filtering criteria. Consider adjusting the thresholds.")
        return None

    filtered.loc[:, 'direction'] = np.where(
        filtered['log2fc'] > 0, "⬆️", "⬇️")
    filtered.loc[:, 'magnitude'] = 2 ** abs(filtered['log2fc'])
    filtered.loc[:, 'expression_change'] = filtered['direction'] + \
        " " + filtered['magnitude'].round(2).astype(str) + "x"

    return filtered
