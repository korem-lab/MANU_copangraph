# Run Episode 2: The Run of the Python
    # Modular script to run the prediction tasks end-to-end (from raw Copan outputs to prediction pickle files)
# Date: 11/04/2024

import glob
import sys
import os
import pandas as pd
import yaml
from PredictionPipelineV3 import *
import argparse

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

def get_X(graph_prefix, out_prefix, nodes, edges, abund=False):

    if abund and not edges:
        print("ABUNDANCE")
        nodes = pd.read_pickle(f'{graph_prefix}.abundance.pkl').T.astype(float)
        nodes_features = nodes.columns.to_series()
        nodes_features.index = range(len(nodes_features))
        nodes_features.to_csv(f'{out_prefix}.abndundance.feature_map.csv', index=None, header=None)
        # Merge nodes and edges into a single dataframe
        X = pd.concat([nodes], axis=1)

    elif abund and edges:
        print("ABUNDANCE+EDGES")
        nodes = pd.read_pickle(f'{graph_prefix}.abundance.pkl').T.astype(float)
        nodes_features = nodes.columns.to_series()
        nodes_features.index = range(len(nodes_features))
        nodes_features.to_csv(f'{out_prefix}.abndundance.feature_map.csv', index=None, header=None)

        edges = pd.read_csv(f'{graph_prefix}.ecolor.csv', index_col=0).T.astype(float)
        edges_features = edges.columns.to_series()
        edges_features.index = range(len(edges_features))
        edges_features.to_csv(f'{out_prefix}.ecolor.feature_map.csv', index=None, header=None)
        assert(all(nodes.index == edges.index))
        X = pd.concat([nodes, edges], axis=1)
    
    elif (nodes and edges):
        print("NODES+EDGES")
        nodes = pd.read_csv(f'{graph_prefix}.ncolor.csv', index_col=0).T.astype(float)
        edges = pd.read_csv(f'{graph_prefix}.ecolor.csv', index_col=0).T.astype(float)
        nodes_features = nodes.columns.to_series()
        edges_features = edges.columns.to_series()
        nodes_features.index = range(len(nodes_features))
        edges_features.index = range(len(edges_features))
        print(">> Number of nodes:", len(nodes.columns))
        nodes_features.to_csv(f'{out_prefix}.ncolor.feature_map.csv', index=None, header=None)
        print(">> Number of edges:", len(edges.columns))
        edges_features.to_csv(f'{out_prefix}.ecolor.feature_map.csv', index=None, header=None)
        # Merge nodes and edges into a single dataframe
        assert(all(nodes.index == edges.index))
        X = pd.concat([nodes, edges], axis=1)

    elif nodes and not edges:
        print("NODES")
        nodes = pd.read_csv(f'{graph_prefix}.ncolor.csv', index_col=0).T.astype(float)
        nodes_features = nodes.columns.to_series()
        nodes_features.index = range(len(nodes_features))
        #nodes.columns = ['n'+str(c) for c in range(len(nodes_features))]
        print(">> Number of nodes:", len(nodes.columns))
        nodes_features.to_csv(f'{out_prefix}.ncolor.feature_map.csv', index=None, header=None)
        # Merge nodes and edges into a single dataframe
        X = pd.concat([nodes], axis=1)

    elif edges and not nodes:
        print("EDGES")
        edges = pd.read_csv(f'{graph_prefix}.ecolor.csv', index_col=0).T.astype(float)
        edges_features = edges.columns.to_series()
        edges_features.index = range(len(edges_features))
        #edges.columns = ['e'+str(c) for c in range(len(edges_features))]
        print(">> Number of edges:", len(edges.columns))
        edges_features.to_csv(f'{out_prefix}.ecolor.feature_map.csv', index=None, header=None)
        # Merge nodes and edges into a single dataframe
        X = pd.concat([edges], axis=1)
    return X


def main(input_dir, graph_name, output_dir, task_name, labels, robust, run_params, model_params, dry_run, nodes, edges, abund, iterations=250, ofld=10, ifld=5, corr=1):
    

    # Load metadata
    labels = get_labels(labels)
    
    # Load features, taken directly from Copan's output
    X = get_X(os.path.join(input_dir, graph_name), os.path.join(output_dir, task_name), nodes, edges, abund)
    X.index = X.index.astype(str)

    # Assert X and md index match (except for the old metadata)
    print(X.index)
    print(labels.index)
    assert(all(X.index == labels.index))
    # Print out the balance of classes
    print(">> Class balance:", labels.outcome.value_counts())

    if remove_correlated_at is not None:
        X = remove_correlated_at(X, corr)

    # Save X and md to output directory as pickle files
    X.to_pickle(f'{output_dir}/{task_name}.X.pkl')
    labels.to_pickle(f'{output_dir}/{task_name}.labels.pkl')

    print(f">> Running with robust={robust}")

    # sort the X vals
    #X = X[X.apply(lambda col: ''.join(col.astype(int).astype(str)), axis=0).sort_values().index]
    
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
        run_name=task_name,
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
        #group_by='group', no groups for moms-pi
        target_col='outcome',
        on_cluster=True,
        seed=73,
        n_jobs=100,
        mem = '4G',
        hours = 24,
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
        conda_env='/burg/pmg/users/az2732/conda_envs/pp3_burg4'
    )
    pred.notes = """
    """

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Copan prediction, with user-defined inputs')
    parser.add_argument('-i', help='Graph output directory')
    parser.add_argument('-g', help='Graph name')
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
    parser.add_argument('--nodes', action='store_true', help='Dry run the script without executing the prediction')
    parser.add_argument('--edges', action='store_true', help='Dry run the script without executing the prediction')
    parser.add_argument('--abund', action='store_true', help='Dry run the script without executing the prediction')
    args = parser.parse_args()
    
    main(
        input_dir=args.i, 
        graph_name=args.g, 
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
        nodes=args.nodes, 
        edges=args.edges, 
        abund=args.abund
    )
