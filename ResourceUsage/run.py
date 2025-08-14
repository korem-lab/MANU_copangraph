import sys
import os
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon, mannwhitneyu
import matplotlib.pyplot as plt
plt.figure(figsize=(3,3), dpi=1400)


DATA = '../data/ResourceUsage'
HRS = 3600.0
GB = 2**20


if __name__ == '__main__':
    
    df = pd.read_csv(os.path.join(DATA, 'coassembly_resource.csv'))
    mwu = pd.DataFrame(columns=['tool_a', 'tool_b', 'statistic', 'pvalue', 'metric', 'coasm_sz'])
    coldict = {'megahit':'red', 'metaspades':'green'}
    asterisk_y = {'megahit': -1, 'metaspades': -2}
    
    # time plot
    plt.clf()
    ax = sns.lineplot(x=df.coasm_sz, y=df['time(s)']/HRS, hue=df.tool, estimator='median')
    sns.scatterplot(x=df.coasm_sz, y=df['time(s)']/HRS, hue=df.tool, ax=ax)
    y_max = df['time(s)'].max() / HRS + 2
    for coasm_sz, g in df.groupby(by=['coasm_sz']):
        tools = set(df.tool) - {'copangraph'}
        for tool in tools:
            copan_vals = g.loc[g.tool == 'copangraph', 'time(s)'].values
            coasm_vals = g.loc[g.tool == tool, 'time(s)'].values
            if len(coasm_vals) == 0:
                continue
            print(coasm_sz, tool)
            print(copan_vals)
            print(coasm_vals)
            res = mannwhitneyu(copan_vals, coasm_vals, alternative='less')
            mwu.loc[len(mwu), :] = ['copangraph', tool, res.statistic, res.pvalue, 'time', coasm_sz[0]]
            print(mwu)
            symbol = ''
            if res.pvalue <= 0.05:
                symbol = '*'
                print(symbol)
            ax.text(int(coasm_sz[0]), y_max + asterisk_y[tool], symbol, color=coldict[tool])
    plt.ylim((0, y_max))
    plt.tight_layout()
    plt.savefig(os.path.join(DATA, 'coassembly_time_labels.pdf'), dpi=1400)
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticklabels([])
    frame1.legend().set_visible(False)
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(DATA, 'coassembly_time.pdf'), dpi=1400)
    
    # mem plot
    asterisk_y = {'megahit': -10, 'metaspades': -20}
    plt.clf()
    ax = sns.lineplot(x=df.coasm_sz, y=df['mem(kb)']/GB, hue=df.tool, estimator='median')
    sns.scatterplot(x=df.coasm_sz, y=df['mem(kb)']/GB, hue=df.tool, ax=ax)
    y_max = df['mem(kb)'].max() / GB + 20
    for coasm_sz, g in df.groupby(by=['coasm_sz']):
        tools = set(df.tool) - {'copangraph'}
        for tool in tools:
            copan_vals = g.loc[g.tool == 'copangraph', 'mem(kb)'].values
            coasm_vals = g.loc[g.tool == tool, 'mem(kb)'].values
            if len(coasm_vals) == 0:
                continue
            print(coasm_sz, tool)
            print(copan_vals)
            print(coasm_vals)
            res = mannwhitneyu(copan_vals, coasm_vals, alternative='less')
            mwu.loc[len(mwu), :] = ['copangraph', tool, res.statistic, res.pvalue, 'mem', coasm_sz[0]]
            print(mwu)
            symbol = ''
            if res.pvalue <= 0.05:
                symbol = '*'
                print(symbol)
            ax.text(int(coasm_sz[0]), y_max + asterisk_y[tool], symbol, color=coldict[tool])
    plt.ylim((0, y_max))
    plt.tight_layout()
    plt.savefig(os.path.join(DATA, 'coassembly_mem_labels.pdf'), dpi=1400)
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticklabels([])
    frame1.legend().set_visible(False)
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(DATA, 'coassembly_mem.pdf'), dpi=1400)
    mwu.to_csv(DATA + 'mwu.csv')
