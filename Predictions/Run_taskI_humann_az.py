# Run Humann on Task I
    # Date: 05/16/2025

import glob
import sys
import os
import pandas as pd
import yaml
from PredictionPipelineV3 import *
import argparse
import re

def get_labels(labels, groups=False):
    df = pd.read_csv(labels, index_col=0)
    df.index = df.index.astype(str)
    assert 'outcome' in df.columns[-1].lower(), "Last column should be outcome"
    if groups:
        assert 'group' in df.columns[-2].lower(), 'Second last column should be group'
    df = df.rename(columns={df.columns[-1]: 'outcome'})  # outcome column can be named differently
    # assert there are no NaN values in the outcome column
    assert df['outcome'].isnull().sum() == 0
    return df

def remove_correlated_at(X, corr=1):
    assert corr == 1
    print('X shape prior to clustering: ', X.shape)
    if corr == 1:
        X = X.loc[:, ~X.T.duplicated()]
        print('X shape post clustering: ', X.shape)
        return X
    # possibly implement clustering at arbitrary correlation


def get_X(labels, organism_specific):

    # Read in humann table
    df = pd.read_csv('/manitou/pmg/projects/korem_lab/Projects/ACU_PLT/mmmbp2/df_path2/TaskI_genefamilies.csv', index_col=0)
    X = df.T
    
    if organism_specific:
        X = X.loc[:, X.columns.str.contains('\|') | (X.columns == 'UNMAPPED')]
    else:
        X = X.loc[:, ~X.columns.str.contains('\|')]
    X.columns = [re.sub(r'[{}[\]":,]', '', col) for col in X.columns]  # remove forbidden characters from column names (for LightGBM)
    
    # X contains every ACU sample, so subselect based on the task-specific labels
    shared_samples = X.index.intersection(labels.index)
    print("Shared samples:", len(shared_samples))
    X = X.loc[shared_samples]
    labels = labels.loc[shared_samples]

    return X, labels

def remove_correlated_at(X, corr=1):
    assert corr == 1
    print('X shape prior to clustering: ', X.shape)
    if corr == 1:
        X = X.loc[:, ~X.T.duplicated()]
        print('X shape post clustering: ', X.shape)
        return X

def main(output_dir, task_name, labels, robust, run_params, model_params, dry_run, iterations=250, ofld=10, ifld=5):
    

    # Load metadata
    labels = get_labels(labels)
    
    # Load features
    ORGANISM_SPECIFIC = False
    X, labels = get_X(labels=labels, organism_specific=ORGANISM_SPECIFIC)

    # Assert X and md index match (except for the old metadata)
    assert(all(X.index == labels.index))
    # Print out the balance of classes
    print(">> Class balance:", labels.outcome.value_counts())
    
    # remove identical features
    X = remove_correlated_at(X, corr=1)

    # Save X and md to output directory as pickle files
    #X.to_pickle(f'{output_dir}/{task_name}.X.pkl')
    #labels.to_pickle(f'{output_dir}/{task_name}.labels.pkl')

    print(f">> Running with robust={robust}")
    
    # === Run Prediction ===
    assert run_params is not None
    run_params = yaml.load(open(run_params), yaml.SafeLoader)
    run_params['feature_selection_method'] = [[FeatureSelectionMethods.LASSO]] if run_params['feature_selection_method'] == 'lasso' else [[FeatureSelectionMethods.NONE]]

    if model_params is None:
        model_params = Defaults.LIGHTGBM_NEW
    else:
        model_params = yaml.load(open(model_params), yaml.SafeLoader)
        print(model_params)

    pred = Prediction(
        X=X,
        md=labels,
        run_name=f'humann.orgspc={ORGANISM_SPECIFIC}.'+task_name,
        model_type=ModelType.LIGHTGBM,
        run_hp_space=run_params,
        model_hp_space=model_params
    )
    print(iterations)

    pred.run(
        inner_folds=ifld,
        outer_folds=ofld,
        robust=robust,
        iterations=iterations,
        stratified=True,
        stratify_by=['outcome'],
        group_by='group',
        target_col='outcome',
        on_cluster=True,
        seed=73,
        n_jobs=150,
        mem = '36G',
        hours = 3,
        cpus=2,
        out_file=f'{output_dir}/humann.orgspc={ORGANISM_SPECIFIC}.{task_name}_1000iters.prediction.pkl',
        save_features=False,
        use_smart_hp_sampler=False,
        sparse_X=True,
        dry_run=dry_run,
        add_name_suffix_to_file=False,
        stagger_seconds=10,
        tryrerun=True,
        run_nested_evaluation=True,
        run_shap=False,
        conda_env='/burg/pmg/users/az2732/conda_envs/pp3_burg4',
        path_to_qp_recover='/burg/pmg/users/ic2465/Projects/MANU_copangraph/Predictions/humann.orgspc=False.TaskI_hmn_p3_runner__fd376fa1c7f748a0a69a28d54ae95f41'
    )
    pred.notes = """
    """

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Copan prediction, with user-defined inputs')
    parser.add_argument('-o', help='Out directory for prediction')
    parser.add_argument('-p', help='Prediction Task name')
    parser.add_argument('--labels', help='Y: the labels and groups')
    parser.add_argument('--robust', type=int, default=3, help='Robust parameter for pp3. Default=1.')
    parser.add_argument('--run_params', default=1, help='run params')
    parser.add_argument('--model_params', default=-1, help='model params')
    parser.add_argument('--iterations', type=int, default=250, help='model params')
    parser.add_argument('--ifld', type=int, default=5, help='model params')
    parser.add_argument('--ofld', type=int, default=10, help='model params')
    parser.add_argument('--dryrun', action='store_true', help='Dry run the script without executing the prediction')
    args = parser.parse_args()
    
    main(
        output_dir=args.o, 
        task_name=args.p, 
        labels=args.labels, 
        robust=args.robust, 
        run_params=args.run_params, 
        model_params=args.model_params,
        iterations=args.iterations,
        ofld=args.ofld,
        ifld=args.ifld,
        dry_run=args.dryrun, 
    )
