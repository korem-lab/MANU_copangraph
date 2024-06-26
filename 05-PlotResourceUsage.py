import sys
import os
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

plt.figure(figsize=(4,4))

DATA = './data/ResourceUsage'
HRS = 3600.0
GB = float(2**20)


if __name__ == '__main__':
    
    df = pd.read_csv(os.path.join(DATA, 'coassembly_resource.csv'))
    wcox = pd.DataFrame(columns=['tool_a', 'tool_b', 'statistic', 'pvalue', 'metric'])
    
    # time plot
    plt.clf()
    ax = sns.lineplot(x=df.coasm_sz, y=df['time(s)']/HRS, hue=df.tool, estimator='median', errorbar=None)
    sns.scatterplot(x=df.coasm_sz, y=df['time(s)']/HRS, hue=df.tool, ax=ax, legend=None)
    y_max = df['time(s)'].max() / HRS + 5
    #for coasm_sz, g in df.groupby(by=['coasm_sz']):
    #    coasm_sz = coasm_sz[0]
    #    copan_vals = g.loc[g.tool == 'copangraph', 'time(s)'].values
    #    coasm_vals = g.loc[g.tool == 'megahit', 'time(s)'].values
    #    res = mannwhitneyu(copan_vals, coasm_vals, alternative='less')
    #    wcox.loc[len(wcox), :] = ['copangraph', 'megahit', res.statistic, res.pvalue, 'time']
    #    symbol = ''
    #    print(res.pvalue)
    #    if res.pvalue <= 0.05:
    #        symbol = '*'
    #    ax.text(coasm_sz, y_max, symbol, color='black')
    plt.tight_layout()
    plt.savefig(os.path.join(DATA, 'coassembly_time.pdf'), dpi=1400, bbox_inches='tight')
    
    # mem plot
    plt.clf()
    y_max = df['time(s)'].max() / GB + 20
    ax = sns.lineplot(x=df.coasm_sz, y=df['mem(kb)']/GB, hue=df.tool, estimator='median', errorbar=None)
    sns.scatterplot(x=df.coasm_sz, y=df['mem(kb)']/GB, hue=df.tool, ax=ax, legend=None)
    #for coasm_sz, g in df.groupby(by=['coasm_sz']):
    #    coasm_sz = coasm_sz[0]
    #    copan_vals = g.loc[g.tool == 'copangraph', 'mem(kb)'].values
    #    coasm_vals = g.loc[g.tool == 'megahit', 'mem(kb)'].values
    #    res = mannwhitneyu(copan_vals, coasm_vals, alternative='less')
    #    wcox.loc[len(wcox), :] = ['copangraph', 'megahit', res.statistic, res.pvalue, 'memory']
    #    symbol = ''
    #    print(res.pvalue)
    #    if res.pvalue <= 0.05:
    #        symbol = '*'
    #    ax.text(coasm_sz, y_max, symbol, color='black')
    plt.tight_layout()
    plt.savefig(os.path.join(DATA, 'coassembly_mem.pdf'), dpi=1400, bbox_inches='tight')