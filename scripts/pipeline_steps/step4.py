import xgboost as xgb
import pandas as pd
from scripts.globals import XGB_MODEL
from scripts.config import FILTER_THRESHOLD


def reorder_columns_for_prediction(df):

    id_array = df.pop("id")

    column_order = ['pred_energy',
                    'pred_num_basepairs',
                    'pred_seed_basepairs',
                    'ta_log10',
                    'sps_mean',
                    'anchor_a',
                    '6mer_seed',
                    'match_8',
                    '6mer_seed_1_mismatch',
                    'compensatory_site',
                    'supplementary_site',
                    'supplementary_site_2',
                    'empty_seed',
                    '9_consecutive_match_anywhere',
                    'mirna_conservation',
                    'seed_8mer',
                    'seed_7mer_a1',
                    'seed_7mer_m8',
                    'seed_compensatory',
                    'seed_clash_2',
                    'seed_clash_3',
                    'seed_clash_4',
                    'seed_clash_5',
                    'mre_au_content',
                    'local_au_content']

    df = df.reindex(columns=column_order)

    return df, id_array


def make_predictions_with_xgb(df):
    model = xgb.Booster()
    model.load_model(XGB_MODEL)

    data_matrix = xgb.DMatrix(df)
    return model.predict(data_matrix)


def create_results_df(id_array, predictions, filter_range=FILTER_THRESHOLD):
    df = pd.DataFrame({'id': id_array, 'prediction': predictions})
    df[['id', 'is_mutated']] = df['id'].str.rsplit('_', n=1, expand=True)
    df["is_mutated"] = df["is_mutated"].eq("mut")

    pivot_df = df.pivot(index='id', columns='is_mutated',
                        values='prediction').reset_index()
    pivot_df.columns = ['id', 'wt_prediction', 'mut_prediction']
    pivot_df["pred_difference"] = (
        pivot_df["mut_prediction"] - pivot_df["wt_prediction"]).round(3)

    # filter out predictions that are outside of the quantile range
    mask = (pivot_df["pred_difference"].abs() >= filter_range)
    return pivot_df[mask]
