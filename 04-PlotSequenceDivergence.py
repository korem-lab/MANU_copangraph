import numpy as np
import pandas as pd
import sys
import os
import glob
import re
import seaborn as sns
from scipy.stats import percentileofscore
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

DATA = './data/SequenceDivergence/'
SEQDIVS = [0.0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2]
figure(figsize=(6,2), dpi=600)

def plot_kde(data, sd, name, min_x, within=True):
    kde = sns.kdeplot(data)
    kde_x, kde_y = kde.get_lines()[0].get_data()
    plt.fill_between(kde_x, kde_y, where=(kde_x >= 1-sd if within else kde_x <= 1-sd), interpolate = True, color='darkblue', alpha=0.50)
    plt.axvline(1.0-sd, color='darkred', linestyle='--')
    plt.axvline(data.mean(), color='darkblue', linestyle='--')
    plt.xlim((min_x, 1.0))
    plt.xlabel(f'{1.0-sd:.3f}% threshold')
    plt.tight_layout()
    plt.savefig(name)
    plt.clf()

    
    
    
if __name__ == '__main__':
    dat = glob.glob(os.path.join(DATA, '*.npy'))
    within_mean = {re.findall('megahit_20_(0.[0-9]+)_perc', e)[0]: np.load(e) for e in dat if ('within' in e and 'mean' in e)}
    within_min = {re.findall('megahit_20_(0.[0-9]+)_perc', e)[0]: np.load(e) for e in dat if ('within' in e and 'min' in e)}
    between = {re.findall('megahit_20_(0.[0-9]+)_perc', e)[0]: np.load(e) for e in dat if 'between' in e}
    
    df = pd.DataFrame(columns=['measure', 'sd', 'percentile_of_score']) 
    for sd in SEQDIVS:
        sd_str = str(sd)
        print(sd_str)
        plot_kde(within_mean[sd_str], sd, os.path.join(DATA, f'within_mean_{sd_str}.pdf'), 0.95, within=True)
        df.loc[len(df), :] = ['within_mean', f'{1-sd:.3f}', 100 - percentileofscore(within_mean[sd_str], 1-sd)]
        plot_kde(within_min[sd_str], sd, os.path.join(DATA, f'within_min_{sd_str}.pdf'), 0.95, within=True)
        df.loc[len(df), :] = ['within_min', f'{1-sd:.3f}', 100 - percentileofscore(within_min[sd_str], 1-sd)]
        plot_kde(between[sd_str], sd, os.path.join(DATA, f'between_{sd_str}.pdf'), 0.5, within=False)
        df.loc[len(df), :] = ['between', f'{1-sd:.3f}', 100 - percentileofscore(between[sd_str], 1-sd)]
    df.to_csv(os.path.join(DATA, 'percentile_of_score.csv'))
        
        