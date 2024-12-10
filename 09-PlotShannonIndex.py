import sys

import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import percentileofscore
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import  figure
import pandas
import numpy
import os
import glob
import re

figure(figsize=(3,3),dpi=1400)

DATA='/manitou/pmg/users/ic2465/Projects/MANU_copangraph_2022/analyses/homology_evaluation/results_before_pe_ext_separation'
OUTDIR='./data/TaxonomicResolution/'

def build_plot_df(files):
    records = [pd.read_pickle(e) for e in files]
    summary_stats = list()
    template = pd.DataFrame(
        index=range(8), columns=['sample_id', 'depth', 'asm', 'dt', 'tax_lvl', 'non_zero_filtered', 'summary', 'value'],
        dtype=object
    )
    for r in records:
        sample_id = int(re.findall('sample_([0-9]+)_', r.sample_id[0])[0])

        # fill record info
        template.index = r.index
        template.sample_id = sample_id
        template.depth = r.depth
        template.asm = r.asm
        template['dt'] = r['dt']
        template.tax_lvl = r.index
        filt = [type(e) == int for e in r.columns]
        data = r.loc[:, filt]

#        # fill sum, mean, median
#        template.summary = 'sum'
#        template.non_zero_filtered = False
#        template.value = data.sum(axis=1)
#        summary_stats.append(template.copy())
#
        template.summary = 'mean'
        template.non_zero_filtered = False
        template.value = data.mean(axis=1)
        summary_stats.append(template.copy())
#
#        template.summary = 'median'
#        template.non_zero_filtered = False
#        template.value = data.median(axis=1)
#        summary_stats.append(template.copy())

        template.summary = 'mean'
        template.non_zero_filtered = True
        template.value = data.apply(lambda x: x[x!=0].mean(), axis=1)
        summary_stats.append(template.copy())

#        template.summary = 'median'
#        template.non_zero_filtered = True
#        template.value = data.apply(lambda x: x[x!=0].median(), axis=1)
#        summary_stats.append(template.copy())

    final_df = pd.concat(summary_stats)
    final_df.index = range(final_df.shape[0])
    return final_df

def plot_figs(plot_df, plot_type, asm='megahit', non_zero_filtered=False):

    # filter for assembler
    plot_df = plot_df.loc[plot_df.asm == asm, :]
    plot_df = plot_df.loc[(plot_df.tax_lvl == 'Species') | (plot_df.tax_lvl == 'Genera') | (plot_df.tax_lvl == 'Family'), :]

    # filter for zero
    plot_df = plot_df.loc[plot_df.non_zero_filtered == non_zero_filtered, :]

    filtered = 'non_zero_filtered' if non_zero_filtered else 'unfiltered'
    # plot mean
    dat = plot_df.loc[plot_df.summary == 'mean', :]
    ax = sns.lineplot(x=dat['dt'], y=dat.value, hue=dat.tax_lvl,
                      estimator='mean', errorbar='sd', legend=True)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, f'{plot_type}{asm}_shannon_index_labels.pdf'))
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticklabels([])
    frame1.legend().set_visible(False)
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR,f'{plot_type}{asm}_shannon_index_nolabels.pdf'))
    plt.clf()

if __name__ == '__main__':

    all_files = glob.glob(os.path.join(DATA, '*.pkl'))
    all_files = [e for e in all_files if '10M' in e]
    plot_type = 'entropies'
    plot_df = build_plot_df([e for e in all_files if plot_type in e])
    plot_figs(plot_df, plot_type, asm='megahit', non_zero_filtered=False)
