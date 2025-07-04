# !/usr/bin/python
# -*- coding: utf-8 -*-
# Author:lihuiru

import os
import sys
from pathlib import Path
import joblib
import pandas as pd
from joblib import load
from pandas.errors import EmptyDataError
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier

MAPPING_DICT = {0: 'avian', 1: 'human'}

def transform_df(df):
    df = df.fillna('')

    # Group by 'Strain ID' and 'Protein Type', then aggregate 'Adaptation Markers'
    transformed = df.groupby(['Strain ID', 'Protein Type']).agg({
        'Adaptation Markers': lambda x: ','.join(x.astype(str))
    }).reset_index()

    # Calculate the count of 'Adaptation Markers'
    transformed['Number of Adaptation Markers'] = transformed['Adaptation Markers'].apply(
        lambda x: len([item for item in x.split(',') if item != '']))

    # Adding empty columns for 'Sequence Type' since it is not provided in the original data
    transformed['Sequence Type'] = ''

    return transformed


def explode_markers(df, input_marker_path, output_directory, prefix):
    """
    Explode the markers from the CSV file into individual mutations and create a crosstab matrix.

    Parameters:
    - df: A DataFrame for exploding.
    - input_marker_path: Path to the CSV file containing marker data.
    - output_directory: Directory where the output files will be saved.
    - prefix: Prefix for the output filenames.

    Returns:
    - DataFrame with the crosstab matrix of markers.
    """
    # Extract the filename without extension
    input_marker_filename = os.path.basename(input_marker_path).rsplit(".", 1)[0]

    # Read the CSV file into a DataFrame
    df = transform_df(df)
    df_original = df.copy()

    # Group and sum the 'Number of Adaptation markers' by 'Strain ID'
    grouped = df.groupby('Strain ID')['Number of Adaptation Markers'].sum()

    # Identify strains with no adaptive markers
    mask = df['Strain ID'].map(grouped) == 0
    no_adaptive_marker_count = mask.sum()
    if no_adaptive_marker_count > 0:
        print(f"There are {no_adaptive_marker_count} strains without any markers.")

    # Explode the 'Adaptation markers' column from a comma-separated string into a list
    df['Adaptation Markers'] = df['Adaptation Markers'].str.split(',')

    # Expand the list into a new DataFrame and merge it back
    df = df.explode('Adaptation Markers')
    df['Protein Type'] = df['Protein Type'].apply(lambda x: 'NA' if x == 'N2' else x)
    df['marker_Protein'] = df['Adaptation Markers'] + '_' + df['Protein Type']

    # Create a new DataFrame with columns as possible 'marker_Protein' values, rows as 'Strain ID'
    df_matrix = pd.crosstab(df['Strain ID'], df['marker_Protein'])
    df_matrix = df_matrix.applymap(lambda x: 1 if x > 0 else 0)
    df_matrix.reset_index(inplace = True)
    df_label = pd.merge(df_matrix, df_original, on = 'Strain ID')

    # Drop duplicates based on 'Strain ID'
    df_label.drop(labels = ['Adaptation Markers', 'Number of Adaptation Markers', 'Protein Type'], axis = 1,
                  inplace = True)

    df_label.drop_duplicates(subset = "Strain ID", keep = "first", inplace = True)

    output_matrix_filename = f'{output_directory}/{prefix}{input_marker_filename}_matrix.csv'

    if '_' in df_label.columns:
        df_label.drop(columns = ['_'],inplace = True)
    df_label.to_csv(output_matrix_filename, index = False)

    return df_label


def ensemble_predict(ensemble_container, X):
    all_preds = []
    all_probas = []

    for model, features in zip(ensemble_container['models'],
                               ensemble_container['feature_sets']):
        X_subset = X[features]
        pred = model.predict(X_subset)
        all_preds.append(pred)
        if hasattr(model, 'predict_proba'):
            proba = model.predict_proba(X_subset)[:, 1]
            all_probas.append(proba)
    final_pred = np.round(np.mean(all_preds, axis = 0)).astype(int)
    final_proba = np.mean(all_probas, axis = 0) if all_probas else None

    return final_pred, final_proba
def get_explode_marker_file(input_marker_path, output_directory = ".", prefix = ""):
    """
       Processes input marker file(s) by exploding the adaptive markers into individual columns
       indicating the presence or absence of each marker in each strain. This function handles
       both individual files and directories containing multiple marker files.

       Parameters:
       - input_marker_path (str/Path): The path to the input marker file or directory.
       - output_directory (str): The directory where the output files will be saved.
       - prefix (str): The prefix to be added to the output filenames.

       Returns:
       - DataFrame: The combined DataFrame of all processed marker files, or the single processed file.
    """
    os.makedirs(output_directory, exist_ok = True)

    if Path(input_marker_path).is_dir():
        all_df = pd.DataFrame()
        for root, _, files in os.walk(input_marker_path):
            for filename in files:
                file_path = os.path.join(root, filename)
                if os.stat(file_path).st_size > 0:
                    try:
                        df = pd.read_csv(file_path)
                        strain_id = filename.split('_markers.csv')[0]
                        if df.empty:
                            # print(df)
                            df.loc[0,'Strain ID'] = strain_id
                        all_df = pd.concat([all_df, df])
                    except EmptyDataError:
                        print(f"Skipping empty data file:{file_path}")
                else:
                    print(f"Skipping empty file:{file_path}")
        if not all_df.empty:
            result_df = explode_markers(all_df, "combined_markers.csv", output_directory, prefix)
        else:
            result_df = pd.DataFrame(columns = ['Strain ID'])
            for root, _, files in os.walk(input_marker_path):
                for idx,filename in enumerate(files):
                    strain_id = filename.split('_markers.csv')[0]
                    result_df.loc[idx,'Strain ID'] = strain_id

    elif Path(input_marker_path).is_file():
        if os.stat(input_marker_path).st_size > 0:
            try:
                df = pd.read_csv(input_marker_path)
                result_df = explode_markers(df, input_marker_path, output_directory, prefix)
            except EmptyDataError:
                print(f"Input file is empty:：{input_marker_path}")
                return None
        else:
            print(f"Input file is empty:{input_marker_path}")
            return None
    else:
        print(f"Error: {input_marker_path} is not a valid file or directory", file = sys.stderr)
        return None
    return result_df

def reorder_rows(df, threshold):
    """
    Reorder rows based on the threshold for 'Human' probability and ensure the 'Human' row is
    placed before the 'Avian' row if its probability exceeds the threshold.

    Parameters:
    - df (pd.DataFrame): DataFrame containing prediction results with columns 'Strain ID', 'Host', and 'Probability'.
    - threshold (float): The probability threshold for classifying a strain as 'Human'.

    Returns:
    - pd.DataFrame: A reordered DataFrame with 'Human' rows (if probability > threshold) placed before 'Avian' rows
      for each 'Strain ID'.
    """
    result = []
    unique_strain_ids = df['Strain ID'].unique()
    for strain_id in unique_strain_ids:
        subset = df[df['Strain ID'] == strain_id]
        human_row = subset[(subset['Host'] == 'Human') & (subset['Probability'] > threshold)]
        avian_row = subset[subset['Host'] == 'Avian']
        if not human_row.empty:
            result.append(human_row.iloc[0])
            result.append(avian_row.iloc[0])
        else:
            result.append(avian_row.iloc[0])
            result.append(subset[subset['Host'] == 'Human'].iloc[0])

    return pd.DataFrame(result)


def load_ensemble_model(model_path):
    ensemble = joblib.load(model_path)
    for model in ensemble['models']:
        if isinstance(model, GradientBoostingClassifier):
            if not hasattr(model, 'n_features_'):
                model.n_features_ = model.n_features
    return ensemble

def predict_new_data(input_marker_path, model_path, threshold, output_directory = ".",
                     prefix = ""):
    """
    Predict class labels for new data using a trained model, threshold, and a set of top features.

    Parameters:
    - input_marker_path (str): Path to the input marker file or directory containing marker files.
    - model_path (str): Path to the saved model file.
    - threshold_path (str): Path to the saved threshold file.
    - top_features_path (str): Path to the saved top features file.
    - output_directory (str): Directory where the output files will be saved.
    - prefix (str): Prefix for the output filenames.

    Returns:
    - DataFrame: A DataFrame with strain IDs and their predicted class labels.
    """
    # Process the input marker file first
    add_prefix = prefix + "_" if prefix else ""
    processed_data = get_explode_marker_file(input_marker_path, output_directory, add_prefix)
    processed_data.set_index("Strain ID", inplace = True)
    # Load the trained model
    loaded_ensemble = load(model_path)

    # Initialize a DataFrame with all required features filled with zeros
    all_required_features = []
    for feature_set in loaded_ensemble['feature_sets']:
        all_required_features.extend(feature_set)
    all_required_features = list(set(all_required_features))  # Remove duplicates

    # Create a complete DataFrame with all required features
    complete_data = pd.DataFrame(0, index = processed_data.index, columns = all_required_features)

    # Fill in the values we have from processed_data
    for col in processed_data.columns:
        if col in all_required_features:
            complete_data[col] = processed_data[col]
    # Now predict using the ensemble
    y_pred, new_data_proba = ensemble_predict(loaded_ensemble, complete_data)
    # Calculate the probability of being 'human' or 'avian' for each prediction
    human_probabilities = [round(prob, 3) for prob in new_data_proba]
    avian_probabilities = [round(1-prob, 3) for prob in new_data_proba]
    # Create DataFrame with prediction results and probabilities
    prediction_results = []
    for i in range(len(new_data_proba)):
        strain_id = complete_data.index[i]
        human_prob = human_probabilities[i]
        avian_prob = avian_probabilities[i]
        # Apply threshold to determine the host
        if human_prob >= threshold:
            host = 'Human'
            prob = human_prob
        else:
            host = 'Avian'
            prob = avian_prob
        prediction_results.append({
            'Strain ID': strain_id,
            'Host': host,
            'Probability': prob
        })

    # Convert to DataFrame
    prediction_results_df = pd.DataFrame(prediction_results).reset_index(drop = True)

    # Save to CSV
    prediction_results_df.to_csv(f"{output_directory}/{add_prefix}prediction.csv", index = False, encoding = 'utf-8-sig')

    return prediction_results_df


if __name__ == "__main__":
    pass

