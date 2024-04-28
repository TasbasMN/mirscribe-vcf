import xgboost as xgb
from scripts.globals import XGB_MODEL, QUANTILE_RANGE
import pandas as pd


# def make_predictions_regressor(df, df_filtered):

#     # Create an empty DMatrix model
#     dmatrix_model = xgb.DMatrix(df_filtered)

#     # Load the pre-trained model
#     model = xgb.Booster()
#     model.load_model(XGB_MODEL)

#     # Make predictions on the df_filtered DataFrame
#     predictions = model.predict(dmatrix_model)

#     # Append the predictions to the "prediction" column in the df DataFrame
#     df["prediction"] = predictions
#     df["binary_prediction"] = (df["prediction"] > 0.5).astype(int)
#     df.sort_values(["id", "is_mutated"], inplace=True)

#     # creating wt and mut dfs
#     wt = df[df.is_mutated == 0].reset_index(drop=True)
#     mut = df[df.is_mutated == 1].reset_index(drop=True)

#     # Calculate the difference between wt and mut predictions
#     wt['pred_difference'] = mut['prediction'] - wt['prediction']
#     wt['pred_difference_binary'] = mut['binary_prediction'] - \
#         wt['binary_prediction']

#     # Merge the difference values back to the original df DataFrame
#     df = df.merge(
#         wt[['id', 'pred_difference', 'pred_difference_binary']], on='id', how='left')

#     return df

# def filter_columns_for_xgb_prediction(df):
#     cols_to_keep = [
#         "pred_energy",
#         "pred_num_basepairs",
#         "pred_seed_basepairs",
#         "ta_log10",
#         "sps_mean",
#         "anchor_a",
#         "6mer_seed",
#         "match_8",
#         "6mer_seed_1_mismatch",
#         "compensatory_site",
#         "supplementary_site",
#         "supplementary_site_2",
#         "empty_seed",
#         "9_consecutive_match_anywhere",
#         "mirna_conservation",
#         "seed_8mer",
#         "seed_7mer_a1",
#         "seed_7mer_m8",
#         "seed_compensatory",
#         "seed_clash_2",
#         "seed_clash_3",
#         "seed_clash_4",
#         "seed_clash_5",
#         "mre_au_content",
#         "local_au_content"
#     ]
#     return df[cols_to_keep]


# def make_predictions(df_with_features):
#     """Make predictions using the XGBoost regressor."""
#     df_filtered = filter_columns_for_xgb_prediction(df_with_features)
#     return make_predictions_regressor(df_with_features, df_filtered)



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


def create_prediction_result_df(id_array, predictions, filter_range=QUANTILE_RANGE, brief_output=False):

    df = pd.DataFrame({'id': id_array, 'prediction': predictions})

    # Convert predictions to binary 
    df["binary_prediction"] = df["prediction"].gt(0.5).astype(int)

    # Split 'id' and create 'is_mutated' 
    df[['id', 'is_mutated']] = df['id'].str.rsplit('_', n=1, expand=True)
    df["is_mutated"] = df["is_mutated"].eq("mut")

    # Calculate pred_differences
    result = df.groupby(['id', 'is_mutated'])['prediction'].first().unstack().diff(axis=1).iloc[:, -1]

    # Merge the result back into the original DataFrame
    df = pd.merge(df, result.reset_index(name='pred_difference'), on='id')
    
    df["pred_difference"] = df["pred_difference"].round(3)
    
    # filter out predictions that are outside of the quantile range
    mask = (df["pred_difference"].abs() >= filter_range)
    df = df[mask]
    
    if brief_output:
        return df[['id', 'pred_difference']].drop_duplicates(subset='id').sort_values('id').reset_index(drop=True)

    return df
        