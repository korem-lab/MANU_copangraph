# Predictions on clinical variables for VRE (Task I) for ACU
# Date: 02/22/2025

import glob
import sys
import os
import pandas as pd
from PredictionPipelineV3 import *
import argparse
from CopanACUPredUtils import load_old_metadata


def get_md_from_input_folder(input_dir, task_name):
    df = pd.read_csv(f'{input_dir}/{task_name}.outcome.csv', index_col=0)
    df.index = df.index.astype(str)
    assert 'outcome' in df.columns[1].lower(), "it is expected that outcome column is third column in csv and it has outcome in name"
    df = df.rename(columns={df.columns[1]: 'outcome'})  # outcome column can be named differently
    df['outcome'] = df['outcome'].replace({'Persistence': 1, 'Clearance': 0})

    # assert there are no NaN values in the outcome column
    assert df['outcome'].isnull().sum() == 0

    # assert only 1 and 0 values in the outcome column
    assert set(df['outcome'].unique()) == {0, 1}

    # assert StudyID is present in columns (to group by)
    assert 'StudyID' in df.columns

    return df


def get_X_from_input_folder(input_dir, task_name, output_dir):
    if task_name == 'oldXoldSeq':
        nodes = pd.read_pickle('input/mdro_pos_05_mo1000_ms100.ncolor.pkl').astype(float)
        nodes.columns = ['n'+str(c) for c in nodes.columns]

        edges = pd.read_pickle('input/mdro_pos_05_mo1000_ms100.ecolor.pkl').astype(float)
        edges.columns = ['e'+str(c) for c in edges.columns]

        return pd.concat([nodes, edges], axis=1)

    nodes = pd.read_csv(f'{input_dir}/{task_name}.ncolor.csv', index_col=0).T.astype(float)
    edges = pd.read_csv(f'{input_dir}/{task_name}.ecolor.csv', index_col=0).T.astype(float)
    assert(all(nodes.index == edges.index))

    nodes_features = nodes.columns.to_series()
    edges_features = edges.columns.to_series()
    nodes_features.index = range(len(nodes_features))
    edges_features.index = range(len(edges_features))
    nodes.columns = ['n'+str(c) for c in range(len(nodes_features))]
    edges.columns = ['e'+str(c) for c in range(len(edges_features))]
    print(">> Number of nodes:", len(nodes.columns))
    nodes_features.to_csv(f'{output_dir}/{task_name}.ncolor.feature_map.csv', index=None, header=None)
    print(">> Number of edges:", len(edges.columns))
    edges_features.to_csv(f'{output_dir}/{task_name}.ecolor.feature_map.csv', index=None, header=None)

    # Merge nodes and edges into a single dataframe
    X = pd.concat([nodes, edges], axis=1)

    return X


def main(input_dir, task_name, robust, output_dir=None, specific_metadata=None, task_code=None, dry_run=False):

    if output_dir is None:
        output_dir = input_dir
        print("Warning: Output directory not specified. Saving results in the input directory.")
        if output_dir.startswith('/manitou'):
            raise ValueError("Writing directly to /manitou is prohibited. Change the output directory.")

    if task_code is not None:
        output_dir = f'{output_dir}/{task_code}'

    os.makedirs(output_dir, exist_ok=True)
    
    print(f"=== RUNNING CLINICAL PREDICTION ON TASK: {task_name} ===")
    print("> Output directory:", output_dir)
    
    # Load metadata
    md = pd.read_csv('/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/plt/TaskI/TaskI.outcome.csv', index_col=0)
    md.index = md.index.astype(str)
    md = md.rename(columns={'Outcome': 'outcome'})
    
    # Load clinical features
    X = pd.read_csv("data/TaskI_clincal_merge.csv", index_col=0)
    # convert index to int
    X.index = X.index.astype(int)
    X.index = X.index.astype(str)
    # Select necessary columns for prediction
    X = X[['Age', 'Sex', 'BMIkgm2', 'Race', 'Ethnicity', 'MELDScore', 'ChildPughScore', 'DaysPostTx', 'course_days_since_last_exposure',
           'PrePeriPost', 'reasonhcv', 'reasonhbv', 'reasonarld', 'reasonnafld', 'reasonhcc', 'reasonother']]
    X['Race'] = X['Race'].replace({98.0: 4})
    X['Ethnicity'] = X['Ethnicity'].replace({-99: 0})
    X['PrePeriPost'] = X['PrePeriPost'].replace({'Pre': 0, 'Peri': 1, 'Post': 2})
    for col in X.columns:
        X[col] = X[col].replace(float('inf'), float('nan'))
    categorical_columns = ['Sex', 'Race', 'Ethnicity', 'PrePeriPost', 'reasonhcv', 'reasonhbv', 'reasonarld', 'reasonnafld', 'reasonhcc', 'reasonother']
    X[categorical_columns] = X[categorical_columns].astype('category')

    # Assert X and md index match (except for the old metadata)
    shared_samples = X.index.intersection(md.index)
    X = X.loc[shared_samples]
    md = md.loc[shared_samples]


    # Print out the balance of classes
    print(">> Class balance:")
    print(md['outcome'].value_counts())

    # Save X and md to output directory as pickle files
    X.to_pickle(f'{output_dir}/{task_name}.X.pkl')
    md.to_pickle(f'{output_dir}/{task_name}.md.pkl')

    print(f">> Running with robust={robust}")
    
    # === Run Prediction ===
    run_params = {
        'steps': [
            ['Imputer'],
        ],
        'imputer_method': ['knn'],
    }

    model_params = Defaults.LIGHTGBM_NEW

    RUN_NAME = task_code if task_code is not None else task_name
    pred = Prediction(
        X=X,
        md=md,
        run_name=RUN_NAME,
        model_type=ModelType.LIGHTGBM,
        run_hp_space=run_params,
        model_hp_space=model_params
    )
    pred.run(
        inner_folds=5,
        outer_folds=10,
        robust=robust,
        iterations=250,
        stratified=True,
        stratify_by=['outcome'],  # MUST BE A LIST, otherwise grouping+stratification will not work, only stratification will
        group_by='StudyID',
        target_col='outcome',
        on_cluster=True,
        seed=73,
        n_jobs=20,
        mem = '14G',
        hours = 5,
        cpus=2,
        out_file=f'{output_dir}/{task_name}.prediction.pkl',
        save_features=False,
        use_smart_hp_sampler=False,
        sparse_X=True,
        dry_run=dry_run,
        add_name_suffix_to_file=False,
        stagger_seconds=10,
        tryrerun=True,
        run_nested_evaluation=True,
        run_shap=False,
        conda_env='pp3_burg4'
    )
    pred.notes = """
    """

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run clinical prediction, with user-defined inputs')
    parser.add_argument('-i', help='Input directory where the Copan outputs and metadata files are located')
    parser.add_argument('-x', help='Prefix for the task name (used for both Copan raw file input, metadata input, and prediction pickle output)')
    parser.add_argument('-o', help='Output directory for the prediction results')
    parser.add_argument('--robust', type=int, default=3, help='Robust parameter for pp3. Default=3.')
    parser.add_argument('--dryrun', action='store_true', help='Dry run the script without executing the prediction')
    args = parser.parse_args()

    task_code = "I"
    
    main(input_dir=args.i, task_name=args.x, robust=args.robust, output_dir=args.o, specific_metadata=None, task_code=task_code, dry_run=args.dryrun)
