import os
import sys
import glob
import pandas as pd
import Utils.GraphQualityPlotting as plot

RESULTS_PATH = "./data/GherigGraphQuality/"
ASM = 'megahit'

if __name__ == '__main__':
    
    ss_quality = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_ss_quality.csv'))
    ss_complexity = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_ss_complexity.csv'))
    ss_nX= glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_ss_nX.csv'))
    co_quality = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_co_results.csv'))
    co_complexity = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_co_complexity.csv'))
    co_nX = glob.glob(os.path.join(RESULTS_PATH, f'*{ASM}*_co_nX.csv'))
    
    # write summary files
    plot.write_summary(ss_quality, f'{ASM}_all_ss_quality.csv')
    plot.write_summary(ss_complexity, f'{ASM}_all_ss_complexity.csv')
    plot.write_summary(ss_nX, f'{ASM}_all_ss_nX_.csv')
    plot.write_summary(co_quality, f'{ASM}_all_co_quality.csv')
    plot.write_summary(co_complexity, f'{ASM}_all_co_complexity.csv')
    plot.write_summary(co_nX, f'{ASM}_all_co_nX.csv')
    
    # Plot single-sample distributions 
    plot.ss_recall(RESULTS_PATH, f'{ASM}_ss_recall', ss_quality)
    plot.ss_complexity(RESULTS_PATH, f'{ASM}_complexity', ss_complexity)
    plot.ss_nX(RESULTS_PATH, f'{ASM}_ss_nX', ss_nX)
    
    # Plot multi-sample distributions
    plot.co_recall(RESULTS_PATH, f'{ASM}_co_recall', co_quality)
    plot.co_complexity(RESULTS_PATH, f'{ASM}_co_complexity', co_complexity)
    plot.co_nX(RESULTS_PATH, f'{ASM}_co_nX', co_nX)
    