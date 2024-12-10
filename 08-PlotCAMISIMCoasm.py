import os
import sys
import re
import glob
import pandas as pd
import Utils.GraphQualityPlotting as plot

RESULTS_PATH = "./data/CAMISIMCoassembly/"

if __name__ == '__main__':
    files = glob.glob(os.path.join(RESULTS_PATH, '*_quality.csv'))
    quality_all = pd.concat([pd.read_csv(e) for e in files], ignore_index=True)
    quality_all.to_csv(os.path.join(RESULTS_PATH, 'all_quality.csv'))
    plot.camisim_coasm_quality(RESULTS_PATH, 'cami_coasm', 'cov_F-score', quality_all)
    plot.camisim_coasm_quality(RESULTS_PATH, 'cami_coasm', 'cov_recall', quality_all)
    plot.camisim_coasm_quality(RESULTS_PATH, 'cami_coasm', 'cov_precision', quality_all)
    plot.camisim_coasm_quality(RESULTS_PATH, 'cami_coasm', 'cnx_F-score', quality_all)
    plot.camisim_coasm_quality(RESULTS_PATH, 'cami_coasm', 'cnx_recall', quality_all)
    plot.camisim_coasm_quality(RESULTS_PATH, 'cami_coasm', 'cnx_precision', quality_all)

    
