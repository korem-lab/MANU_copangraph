import os
import sys
import glob
import pandas as pd
from scipy.stats import wilcoxon
from itertools import product
import Utils.GraphQualityPlotting as plot

RESULTS_PATH = "./data/CAMISIMGraphQuality/"
ASM = 'megahit'

if __name__ == '__main__':
    
    quality_df = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_quality.csv'))
    complexity_df = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_complexity.csv'))
    nX_df = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_nX.csv'))

    # concat data 
    quality_df = pd.concat([pd.read_csv(f) for f in quality_df]) 
    complexity_df = pd.concat([pd.read_csv(f) for f in complexity_df]) 
    nX_df = pd.concat([pd.read_csv(f) for f in nX_df]) 
    
    # write summary files
    quality_df.to_csv(os.path.join(RESULTS_PATH, f'{ASM}_all_quality.csv'))
    complexity_df.to_csv(os.path.join(RESULTS_PATH, f'{ASM}_all_complexity.csv'))
    nX_df.to_csv(os.path.join(RESULTS_PATH, f'{ASM}_all_nX_.csv'))
    
    # Plot graph quality metircs by depth
    plot.graph_quality_by_depth(RESULTS_PATH, f'{ASM}_cnx_F-score', quality_df, metric='cnx_F-score')
    plot.graph_quality_by_depth(RESULTS_PATH, f'{ASM}_cov_F-score', quality_df, metric='cov_F-score')
    plot.graph_quality_by_depth(RESULTS_PATH, f'{ASM}_cnx_precision', quality_df, metric='cnx_precision')
    plot.graph_quality_by_depth(RESULTS_PATH, f'{ASM}_cov_precision', quality_df, metric='cov_precision')
    plot.graph_quality_by_depth(RESULTS_PATH, f'{ASM}_cnx_recall', quality_df, metric='cnx_recall')
    plot.graph_quality_by_depth(RESULTS_PATH, f'{ASM}_cov_recall', quality_df, metric='cov_recall')

    # compute Mann-Whitney U
    scores = pd.read_csv(os.path.join(RESULTS_PATH, 'megahit_cnx_F-score_scores.csv'))
    groups = {(d, a):frame for ((d, a), frame) in scores.groupby(by=['depth', 'assembler'])}
    depths = set(scores.depth)
    tools = set(scores.assembler)
    wcox_tests = pd.DataFrame(columns=['tool_a', 'tool_b', 'depth', 'pval', 'statistic', 'hyp'])
    for depth, tool_a, tool_b in product(depths, tools, tools):
        if tool_a != 'copangraph' or tool_a == tool_b:
            continue
        t1 = groups[(depth, tool_a)]
        t2 = groups[(depth, tool_b)]
        assert(all(t1.dataset.values == t2.dataset.values))
        res = wilcoxon(t1.value, t2.value, alternative='greater')
        wcox_tests.loc[len(wcox_tests), :] = [tool_a, tool_b, depth, res.pvalue, res.statistic, 'greater']
    wcox_tests.to_csv(os.path.join(RESULTS_PATH, 'megahit_cnx_F-score_wcox.csv'))
    
    scores = pd.read_csv(os.path.join(RESULTS_PATH, 'megahit_cov_F-score_scores.csv'))
    groups = {(d, a):frame for ((d, a), frame) in scores.groupby(by=['depth', 'assembler'])}
    depths = set(scores.depth)
    tools = set(scores.assembler)
    wcox_tests = pd.DataFrame(columns=['tool_a', 'tool_b', 'depth', 'pval', 'statistic', 'hyp'])
    for depth, tool_a, tool_b in product(depths, tools, tools):
        if tool_a != 'copangraph' or tool_a == tool_b:
            continue
        t1 = groups[(depth, tool_a)]
        t2 = groups[(depth, tool_b)]
        assert(all(t1.dataset.values == t2.dataset.values))
        res = wilcoxon(t1.value, t1.value, alternative='greater')
        wcox_tests.loc[len(wcox_tests), :] = [tool_a, tool_b, depth, res.pvalue, res.statistic, 'greater']
    wcox_tests.to_csv(os.path.join(RESULTS_PATH, 'megahit_cov_F-score_wcox.csv'))
    
    # plot complexity
    #plot.graph_complexity_by_depth(RESULTS_PATH, f'{ASM}_graph_complexity', complexity_df)

    # plot NX
    #plot.graph_NX_by_depth(RESULTS_PATH, f'{ASM}_graph_complexity', nX_df)