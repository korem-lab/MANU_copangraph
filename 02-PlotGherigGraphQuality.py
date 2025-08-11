import os
import sys
import re
import glob
import pandas as pd
import Utils.GraphQualityPlotting as plot

RESULTS_PATH = "/manitou/pmg/projects/korem_lab/Projects/MANU_copangraph/GherigGraphQuality/coassemblies"
#RESULTS_PATH = "/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/mmmbp_all_edgesync"
RESULTS_PATH = "/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/corrected_aag"
#RESULTS_PATH_SS = "/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/subgraph_analysis/MANUSCRIPT_VERSION/"
RESULTS_PATH_SS = "/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/subgraph_analysis/corrected_aag"


CO = 'coasm_cocpg_alllr'
CO = 'coassembly_analysis'
SS = 'ssasm_sscpg_sslr'

if __name__ == '__main__':
    
    #ss_quality = glob.glob(os.path.join(RESULTS_PATH_SS, f'*{SS}*_quality.csv'))
    #ss_complexity = glob.glob(os.path.join(RESULTS_PATH, f'*{SS}*_complexity.csv'))
    #ss_nX= glob.glob(os.path.join(RESULTS_PATH, f'*{SS}*_nX.csv'))
    co_quality = glob.glob(os.path.join(RESULTS_PATH, f'*{CO}*_quality.csv'))
    #print(co_quality)
    #co_complexity = glob.glob(os.path.join(RESULTS_PATH, f'*{CO}*_complexity.csv'))
    #co_nX = glob.glob(os.path.join(RESULTS_PATH, f'*{CO}*_nX.csv'))
    #print(co_quality)
    
    # concat data
    #ss_quality = pd.concat([pd.read_csv(e) for e in ss_quality], ignore_index=True)
    #ss_complexity = pd.concat([pd.read_csv(e) for e in ss_complexity], ignore_index=True)
    #ss_nX = pd.concat([pd.read_csv(e) for e in ss_nX], ignore_index=True)
    co_quality = pd.concat([pd.read_csv(e) for e in co_quality], ignore_index=True)
    #co_complexity = pd.concat([pd.read_csv(e) for e in co_complexity], ignore_index=True)
    #co_nX = pd.concat([pd.read_csv(e) for e in co_nX], ignore_index=True)
   
    # write 
    #ss_quality.to_csv(os.path.join(RESULTS_PATH_SS, f'{SS}_quality_ALL.csv'))
    #ss_complexity.to_csv(os.path.join(RESULTS_PATH, f'{SS}_all_complexity.csv'))
    #ss_nX.to_csv(os.path.join(RESULTS_PATH, f'{SS}_all_nX.csv'))
    co_quality.to_csv(os.path.join(RESULTS_PATH, f'{CO}_quality_ALL.csv'))
    #co_complexity.to_csv(os.path.join(RESULTS_PATH, f'{CO}_all_complexity.csv'))
    #co_nX.to_csv(os.path.join(RESULTS_PATH, f'{CO}_all_nX.csv'))
    
    # Plot single-sample distributions 
    # Recall boxplots across the 10 datasets for coverage connectivity
    #plot.ss_quality(RESULTS_PATH_SS, SS, 'cov_F-score', ss_quality, True)
    #plot.ss_quality(RESULTS_PATH_SS, SS, 'cov_recall', ss_quality, True)
    #plot.ss_quality(RESULTS_PATH_SS, SS, 'cov_precision', ss_quality, True)
    #plot.ss_quality(RESULTS_PATH_SS, SS, 'cnx_F-score', ss_quality, True)
    #plot.ss_quality(RESULTS_PATH_SS, SS, 'cnx_recall', ss_quality, True)
    #plot.ss_quality(RESULTS_PATH_SS, SS, 'cnx_precision', ss_quality, True)
    # Nodes and edges boxplots across the 10 datasets
    #plot.ss_complexity(RESULTS_PATH, f'{ASM}_ss_complexity', ss_complexity)
    # N50 / N90 boxplots
    #plot.ss_nX(RESULTS_PATH, f'{ASM}_ss_nX', ss_nX)
    
    # Plot multi-sample curves
    #  Recall boxplots for 1, 3, 5, 10-sample co-assemblies 
    co_quality.loc[:, 'coasm_sz'] = co_quality.dataset.apply(lambda x: int(re.findall('([0-9]+)_sample', x)[0]))
    plot.co_quality(RESULTS_PATH, CO, 'cov_F-score', co_quality, True)
    plot.co_quality(RESULTS_PATH, CO, 'cov_recall', co_quality, True)
    plot.co_quality(RESULTS_PATH, CO, 'cov_precision', co_quality, True)
    plot.co_quality(RESULTS_PATH, CO, 'cnx_F-score', co_quality, True)
    plot.co_quality(RESULTS_PATH, CO, 'cnx_recall', co_quality, True)
    plot.co_quality(RESULTS_PATH, CO, 'cnx_precision', co_quality, True)

    # macro
    #plot.co_quality(RESULTS_PATH, CO, 'cov_F-score', co_quality)
    #plot.co_quality(RESULTS_PATH, CO, 'cov_recall', co_quality)
    #plot.co_quality(RESULTS_PATH, CO, 'cov_precision', co_quality)
    #plot.co_quality(RESULTS_PATH, CO, 'cnx_F-score', co_quality)
    #plot.co_quality(RESULTS_PATH, CO, 'cnx_recall', co_quality)
    #plot.co_quality(RESULTS_PATH, CO, 'cnx_precision', co_quality)
    #  Node and edge counts for 1, 3, 5, 10-sample co-assemblies 
    #plot.co_complexity(RESULTS_PATH, f'{ASM}_co_complexity', co_complexity)
    #  N50 / N90 counts for 1, 3, 5, 10-sample co-assemblies 
    #plot.co_nX(RESULTS_PATH, f'{ASM}_co_nX', co_nX)
    
