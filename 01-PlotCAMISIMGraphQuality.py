import os
import sys
import glob
import pandas as pd
import Utils.GraphQualityPlotting as plot

RESULTS_PATH = "./data/CAMISIMGraphQuality/"
ASM = 'megahit'

if __name__ == '__main__':
    
    quality_df = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_quality.csv'))
    complexity_df = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_complexity.csv'))
    nX_df = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_nX.csv'))
    
    # write summary files
    plot.write_summary(quality_df, f'{ASM}_all_quality.csv')
    plot.write_summary(complexity_df, f'{ASM}_all_complexity.csv')
    plot.write_summary(nX_df, f'{ASM}_all_nX_.csv')
    
    # Plot graph quality metircs by depth
    plot.graph_quality_plot_coasm(RESULTS_PATH, f'{ASM}_cnx_F-score', quality_df)
    plot.graph_quality_plot_coasm(RESULTS_PATH, f'{ASM}_cov_F-score', quality_df)
    plot.graph_quality_plot_coasm(RESULTS_PATH, f'{ASM}_cnx_precision', quality_df)
    plot.graph_quality_plot_coasm(RESULTS_PATH, f'{ASM}_cov_precision', quality_df)
    plot.graph_quality_plot_coasm(RESULTS_PATH, f'{ASM}_cnx_recall', quality_df)
    plot.graph_quality_plot_coasm(RESULTS_PATH, f'{ASM}_cov_recall', quality_df)
    
    # plot complexity
    plot.graph_complexity_by_depth(RESULTS_PATH, f'{ASM}_graph_complexity', complexity_df)

    # plot NX
    plot.graph_NX_by_depth(RESULTS_PATH, f'{ASM}_graph_complexity', nX_df)