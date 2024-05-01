import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import figures as fc

ASM = 'megahit_20M'
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("usage: <exe> <results_dir>")
        sys.exit()
    
    csvs = glob.glob(os.path.join(sys.argv[1], f'*mo100_ms75_results.csv'))
    dfs = [pd.read_csv(c) for c in csvs]
    df = pd.concat(dfs)
    df.to_csv(os.path.join(sys.argv[1], f'{ASM}_all_results.csv'))
    df.index = range(df.shape[0])
    #fc.resource_plot(df, 'rss')
    #fc.resource_plot(df, 'time')
    fc.graph_quality_plot_depth(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_glocal95_eps40'), 'cnx_F-score', df, filter_on='copangraph', hue_val='dataset')
    #fc.graph_quality_plot_depth(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_glocal95_eps40'), 'cov_F-score', df, filter_on='copangraph', hue_val='dataset')
    #fc.graph_quality_plot_depth(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_glocal95_eps40'), 'cnx_precision', df, filter_on='copangraph', hue_val='dataset')
    #fc.graph_quality_plot_depth(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_glocal95_eps40'), 'cov_precision', df, filter_on='copangraph', hue_val='dataset')
    #fc.graph_quality_plot_depth(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_glocal95_eps40'), 'cnx_recall', df, filter_on='copangraph', hue_val='dataset')
    #fc.graph_quality_plot_depth(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_glocal95_eps40'), 'cov_recall', df, filter_on='copangraph', hue_val='dataset')
    fc.graph_quality_plot_coasm(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_coasm'), 'cnx_F-score', df)
    fc.graph_quality_plot_coasm(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_coasm'), 'cov_F-score', df)
    fc.graph_quality_plot_coasm(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_coasm'), 'cnx_precision', df)
    fc.graph_quality_plot_coasm(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_coasm'), 'cov_precision', df)
    fc.graph_quality_plot_coasm(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_coasm'), 'cnx_recall', df)
    fc.graph_quality_plot_coasm(os.path.join(sys.argv[1], f'pnp-camisim_{ASM}_coasm'), 'cov_recall', df)
    #scores = fc.construct_distribution_of_single_metric_across_datasets('cnx_recall', df)
    #fc.plot_distribution_of_single_metric_across_datasets(scores, 'Recall', 'cnx_recall.pdf')
    #scores.to_csv('cnx_recall.csv')
    #scores = fc.construct_distribution_of_single_metric_across_datasets('cov_recall', df)
    #fc.plot_distribution_of_single_metric_across_datasets(scores, 'Recall', 'cov_recall.pdf')
    #scores.to_csv('cov_recall.csv')

    #scores = fc.construct_distribution_of_single_metric_across_datasets('cnx_precision', df)
    #fc.plot_distribution_of_single_metric_across_datasets(scores, 'Precision', 'cnx_precision.pdf')
    #scores.to_csv('cnx_precision.csv')

    #scores = fc.construct_distribution_of_single_metric_across_datasets('cnx_F-score', df)
    #fc.plot_distribution_of_single_metric_across_datasets(scores, 'F-score', 'cnx_F-score.pdf')
    #scores.to_csv('cnx_F-score.csv')

