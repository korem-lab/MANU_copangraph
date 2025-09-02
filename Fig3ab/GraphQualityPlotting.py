import random
import sys
import os
import re
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
from itertools import combinations, product

#import evaluation.utils.constants as constants
figure(figsize=(3, 3))

ABD_THRESH = 0.05

def graph_quality_by_depth(out_dir, fname, quality_df, metric=None, filter_on=None, hue_val='assembler', base=2, stats=None):

    # ensure inputs are valid
    plt.clf()
    assert metric in [
        'cnx_precision', 'cnx_recall', 'cnx_F-score',
        'cov_precision', 'cov_recall', 'cov_F-score',
    ]

    groupings = ['depth', 'assembler', 'dataset']
    columns = ['depth', 'assembler', 'dataset', 'value']
    
    # group the table
    groups = quality_df.groupby(groupings)
    
    # For each group, compute the relevant metric
    scores = pd.DataFrame(index=range(len(groups)), columns=columns)
    for i, (fields, g) in enumerate(groups):
        scores.iloc[i, :] = (*fields, compute_metric(metric, g))
    scores.to_csv(os.path.join(out_dir, f'{fname}_scores.csv'), index=None)
    groups = {(d, a):frame for ((d, a), frame) in scores.groupby(by=['depth', 'assembler'])}
    depths = set(scores.depth)
    tools = set(scores.assembler)
    stats = pd.DataFrame(columns=['tool_a', 'tool_b', 'depth', 'pval', 'statistic', 'hyp'])
    for depth, tool_a, tool_b in product(depths, tools, tools):
        if tool_a != 'copangraph' or tool_a == tool_b:
            continue
        t1 = groups[(depth, tool_a)]
        t2 = groups[(depth, tool_b)]
        assert(all(t1.dataset.values == t2.dataset.values))
        try:
            res = wilcoxon(t1.value.astype(float).values, t2.value.astype(float).values, alternative='greater')
            stats.loc[len(stats), :] = [tool_a, tool_b, depth, res.pvalue, res.statistic, 'greater']
            res = wilcoxon(t1.value.astype(float).values, t2.value.astype(float).values, alternative='two-sided')
            stats.loc[len(stats), :] = [tool_a, tool_b, depth, res.pvalue, res.statistic, 'two-sided']
        except ValueError:
            stats.loc[len(stats), :] = [tool_a, tool_b, depth, 1, 0, 'greater']
    stats.to_csv(os.path.join(out_dir, f'megahit_{metric}_wcox.csv'))
    # Plot line graph
    if filter_on:
        scores = scores.loc[scores.assembler == filter_on,:]
    ax = sns.lineplot(x=(scores['depth'].astype(np.float64)), y=scores['value'], hue=scores[hue_val],legend=True, errorbar='sd')

    if stats is not None:
        #y_max = 0.5 if 'cnx' in fname else 0.7
        y_pos = scores.groupby(['assembler', 'depth']).apply(lambda x: x.value.mean() + x.value.std()).max()
        coldict = {'megahit_contigs':'green', 'megahit_graph':'red', 'mcvl':'orange'}
        space_dict = {'megahit_contigs':0.00, 'megahit_graph':0.025, 'mcvl':0.05}
        stats = stats.loc[stats.hyp == 'greater',:]
        for i in stats.index:
            tool = stats.loc[i, 'tool_b']
            pval = stats.loc[i, 'pval']
            if pval >= 0.05:
                symbol = ''
            else:
                symbol = '*'
            depth = int(stats.loc[i, 'depth'])
            ax.text(depth, y_pos + space_dict[tool], symbol, color=coldict[tool])
        ax.set_ylim(top=y_pos+0.1)

    # put ticks at exact measurement position 
    ax.set_xlabel('Read pairs')
    ax.set_xscale('log', base=base)
    plt.xticks(ticks=[1*10**6, 2*10**6, 5*10**6, 10*10**6, 20*10**6, 40*10**6, 60*10**6, 80*10**6], labels=[1*10**6, 2*10**6, 5*10**6, 10*10**6, 20*10**6, 40*10**6, 60*10**6, 80*10**6])
    ax.set_ylabel(metric)
    ax.set_title(f'{fname} {metric}')
    #sns.move_legend(ax, 'upper left', bbox_to_anchor=(1, 1))
    fname = os.path.join(out_dir, fname)
    plt.tight_layout()
    plt.savefig(f'{fname}_labels.pdf', dpi=1400, bbox_inches='tight')

    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticklabels([])
    frame1.legend().set_visible(False)
    plt.title('')
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(f'{fname}_nolabels.pdf', dpi=1400, bbox_inches='tight')
    plt.clf()
    
def graph_complexity_by_depth(out_dir, fname, complexity_df):
    """
    Assumes table is:
    depth, dataset, assembler, metric, value
    """
    plt.clf()
    nodes = complexity_df.loc[complexity_df.metric == 'nodes',:]
    nplot = sns.lineplot(x=nodes['depth'].astype(np.float64), y=nodes['value'], hue=nodes['assembler'], legend=True)
    nplot.set_xlabel('Read pairs')
    nplot.set_ylabel('num_nodes')
    nplot.set_title(f'{fname} - nodes')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_nodes.pdf'), bbox_inches='tight')
    plt.clf()

    edges = complexity_df.loc[complexity_df.metric == 'edges',:]
    eplot = sns.lineplot(x=edges['depth'].astype(np.float64), y=edges['value'], hue=edges['assembler'], legend=True)
    eplot.set_xlabel('Read pairs')
    eplot.set_ylabel('num_edges')
    eplot.set_title(f'{fname} - edges')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_edges.pdf'), bbox_inches='tight')
    plt.clf()
    
    
    
def graph_NX_by_depth(out_dir, fname, nX_df):
    """
    Assumes table is:
    depth, dataset, assembler, metric, value
    """
    plt.clf()
    n50 = nX_df.loc[nX_df.metric == 'N50',:]
    n50plot = sns.lineplot(x=n50['depth'].astype(np.float64), y=n50['value'], hue=n50['assembler'], legend=True)
    n50plot.set_xlabel('Read pairs')
    n50plot.set_ylabel('N50')
    n50plot.set_title(f'{fname} - N50')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_N50.pdf'), bbox_inches='tight')
    plt.clf()
    
    n90 = nX_df.loc[nX_df.metric == 'N90',:]
    n90plot = sns.lineplot(x=n90['depth'].astype(np.float64), y=n90['value'], hue=n90['assembler'], legend=True)
    n90plot.set_xlabel('Read pairs')
    n90plot.set_ylabel('N90')
    n90plot.set_title(f'{fname} - N90')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_N90.pdf'), bbox_inches='tight')
    plt.clf()

def ss_quality(out_dir, asm, metric, ss_quality, include_unmapped=True):
    """Assumes table is: 
    dataset, assembler, metric, value
    """
    plt.clf()
    assert metric in [
        'cnx_precision', 'cnx_recall', 'cnx_F-score',
        'cov_precision', 'cov_recall', 'cov_F-score',
    ]
    groups = ss_quality.groupby(['assembler', 'dataset']) 
    # For each group, compute the relevant metric
    scores = pd.DataFrame(index=range(len(groups)), columns=['assembler', 'dataset', 'metric', 'value'])
    for i, (fields, g) in enumerate(groups):
        scores.iloc[i, :] = (*fields, metric, compute_metric(metric, g, include_unmapped))
    scores.loc[:, 'coasm_sz'] = scores.dataset.apply(lambda x: int(re.findall('([0-9])_sample', x)[0]))
    scores.to_csv(os.path.join(out_dir, f'{asm}_{metric}_ss_scores.csv'))
    #sns.boxplot(
    #    showmeans=True,
    #    meanline=True,
    #    meanprops={'color': 'k', 'ls': '-', 'lw': 1.5},
    #    medianprops={'visible': False},
    #    whiskerprops={'visible': False},
    #    zorder=10,
    #    x="quantile_upper_bound",
    #    y="value",
    #    hue='assembler',
    #    data=df,
    #    showfliers=False,
    #    showbox=False,
    #    showcaps=False,
    #    ax=ax
    #)
    scores = scores.loc[scores.coasm_sz == 9, :]
    sns.boxplot(x=scores.assembler, y=scores['value'], showfliers=False, order=['copangraph', 'metacarvel', 'megahit_contigs', 'megahit_graph'], palette='tab10')
    sns.stripplot(x=scores.assembler, y=scores['value'], legend=False, dodge=True, color='black', order=['copangraph', 'metacarvel', 'megahit_contigs', 'megahit_graph'])
    plt.title(f'{metric}', fontsize=10)
    plt.ylim((None, scores['value'].max() + 0.2))
    plt.savefig(os.path.join(out_dir, f'{asm}_ss_{metric}_labels.pdf'), bbox_inches='tight', dpi=1400)

    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticklabels([])
    frame1.legend().set_visible(False)
    plt.title('')
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{asm}_ss_{metric}_nolabels.pdf'), bbox_inches='tight', dpi=1400)
    plt.clf()
    

def camisim_coasm_quality(out_dir, asm, metric, ss_quality):
    """Assumes table is: 
    dataset, assembler, metric, value
    """
    plt.clf()
    assert metric in [
        'cnx_precision', 'cnx_recall', 'cnx_F-score',
        'cov_precision', 'cov_recall', 'cov_F-score',
    ]
    groups = ss_quality.groupby(['assembler', 'dataset']) 
    # For each group, compute the relevant metric
    scores = pd.DataFrame(index=range(len(groups)), columns=['assembler', 'dataset', 'metric', 'value'])
    for i, (fields, g) in enumerate(groups):
        scores.iloc[i, :] = (*fields, metric, compute_metric(metric, g))
    scores.to_csv(f'{asm}_{metric}_scores.csv')
    #sns.boxplot(
    #    showmeans=True,
    #    meanline=True,
    #    meanprops={'color': 'k', 'ls': '-', 'lw': 1.5},
    #    medianprops={'visible': False},
    #    whiskerprops={'visible': False},
    #    zorder=10,
    #    x="quantile_upper_bound",
    #    y="value",
    #    hue='assembler',
    #    data=df,
    #    showfliers=False,
    #    showbox=False,
    #    showcaps=False,
    #    ax=ax
    #)
    sns.boxplot(x=scores.assembler, y=scores['value'], showfliers=False)
    sns.stripplot(x=scores.assembler, y=scores['value'], legend=False, dodge=True, color='black')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{asm}_{metric}.pdf'), bbox_inches='tight')
    plt.clf()
   
def ss_complexity(out_dir, fname, ss_complexity):
    """Assumes table is: 
    dataset, assembler, metric, value
    """
    plt.clf()
    sns.boxplot(
        x=ss_complexity['metric'], y=ss_complexity['value'], hue=ss_complexity['assembler'],
        showFliers=False, legend=True
    )
    sns.stripplot(
        x=ss_complexity['metric'], y=ss_complexity['value'], hue=ss_complexity['assembler'],
        color='black', dodge=True
    )
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_ss_complexity.pdf'), bbox_inches='tight')
    plt.clf()
    
def ss_nX(out_dir, fname, ss_nX):
    """Assumes table is: 
    dataset, assembler, metric, value
    """
    plt.clf()
    sns.boxplot(
        x=ss_nX['metric'], y=ss_nX['value'], hue=ss_nX['assembler'],
        showFliers=False, legend=True
    )
    sns.stripplot(
        x=ss_complexity['metric'], y=ss_nX['value'], hue=ss_nX['assembler'],
        color='black', dodge=True
    )
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_ss_NX.pdf'), bbox_inches='tight')
    plt.clf()
    
def co_quality(out_dir, asm, metric, co_quality, include_unmapped=True):
    """Assumes table is: 
    coassembly, assembler, metric, value
    """
    plt.clf()
    print(co_quality)
    groups = co_quality.groupby(['assembler', 'coasm_sz'])
    assert metric in [
        'cnx_precision', 'cnx_recall', 'cnx_F-score',
        'cov_precision', 'cov_recall', 'cov_F-score',
        'cnx_precision_macro', 'cnx_recall_macro', 'cnx_F-score_macro',
        'cov_precision_macro', 'cov_recall_macro', 'cov_F-score_macro',
    ]
    
    scores = pd.DataFrame(index=range(len(groups)), columns=['assembler', 'coasm_sz', 'metric', 'value'])
    for i, (fields, g) in enumerate(groups):
        scores.iloc[i, :] = (*fields, metric, compute_metric(metric, g, include_unmapped))
    scores.to_csv(os.path.join(out_dir, f'{asm}_{metric}_co_scores.csv'))
    ax=sns.lineplot(x=(scores['coasm_sz'].astype(np.float64)), y=scores['value'], hue=scores['assembler'], legend=True)
    plt.xticks(ticks=[3,6,9], labels=[3,6,9])
    plt.title(f'{metric}', fontsize=10)
    if 'precision' in metric:
        ax.set_ylim(bottom=0.75)
    plt.savefig(os.path.join(out_dir, f'{asm}_co_{metric}_labels.pdf'), dpi=1400, bbox_inches='tight')
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticklabels([])
    frame1.legend().set_visible(False)
    plt.title('')
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{asm}_co_{metric}_nolabels.pdf'), bbox_inches='tight', dpi=1400)
    plt.clf()

def co_complexity(out_dir, fname, co_complexity):
    """Assumes table is: 
    coassembly, assembler, metric, value
    """
    plt.clf()
    nodes = co_complexity.loc[co_complexity.metric == 'nodes',:]
    sns.lineplot(x=nodes['coassembly'].astype(np.float64), y=nodes['value'], hue='assembler')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_nodes.pdf'), dpi=1400, bbox_inches='tight')
    plt.clf()
    edges = co_complexity.loc[co_complexity.metric == 'edges',:]
    sns.lineplot(x=edges['coassembly'].astype(np.float64), y=edges['value'], hue='assembler')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_edges.pdf'), dpi=1400, bbox_inches='tight')
    plt.clf()
    
def co_nX(out_dir, fname, co_nX):
    plt.clf()
    n50 = co_nX.loc[co_nX.metric == 'N50',:]
    sns.lineplot(x=n50['coassembly'].astype(np.float64), y=n50['value'], hue='assembler')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_N50.pdf'), dpi=1400, bbox_inches='tight')
    plt.clf()
    
    n90 = co_nX.loc[co_nX.metric == 'N90',:]
    sns.lineplot(x=n90['coassembly'].astype(np.float64), y=n90['value'], hue='assembler')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'{fname}_N90.pdf'), dpi=1400, bbox_inches='tight')
    plt.clf()


def complexity_plot(fpath, depth, tool, cat):
    plt.clf()
    data = pd.read_csv(fpath)
    tools = {t for t in data.assembler if 'COMET' in t}
    tools.add(f'{tool}-graph')
    tools = sorted(list({t for t in tools if tool in t}))
    filt = data.assembler.isin(tools) & (data.depth == depth)
    data = data.loc[filt, :]
    sns.boxplot(data=data, x=data.assembler, y=data[cat], order=tools)
    ax = sns.stripplot(data=data, x=data.assembler, y=data[cat], order=tools, color='black', jitter=True)
    ax.set_ylabel(f'{cat}')
    ax.set_title(f'{cat} at increasing homology with {depth} reads')
    ax.set(xticklabels=[])
    plt.xticks(rotation=45)
    #plt.tight_layout()
    plt.savefig(f'./{cat}_{depth}_{tool}_complexity.pdf', dpi=1400, bbox_inches='tight')
    plt.show()


def resource_plot(data, resource):
    plt.clf()
    data = data.sort_values(by=['num_samples', 'tool'])
    print(data)
    _resource = 'rss (GB)' if resource == 'rss' else 'time (min)'
    ax = sns.lineplot(x=data.num_samples.astype(int), y=data[_resource], hue=data.tool.astype(str), estimator='median', legend=True, ci=None)
    sns.scatterplot(x=data.num_samples.astype(int), y=data[_resource], hue=data.tool.astype(str), legend=False, ax=ax)
    plt.ylabel(resource)
    plt.xlabel('N samples')
    plt.tight_layout()
    plt.savefig(f'./{resource}_benchmark.pdf', dpi=1400, bbox_inches='tight')
    plt.show()


#def construct_evaluation_dataframe(dataset, analysis, platform, prefix=''):
#    """
#    Concatenates results dataframes from all relevant evaluations
#    """
#
#    evals = list()
#    # construct the keys
#    for ds in constants.DATASETS:
#
#        # Skip over irrelevant datasets
#        if not ds.startswith(dataset):
#            continue
#
#        for depth in constants.DEPTHS:
#            key = constants.make_key(ds, depth)
#            eval_path = constants.analysis_path(key, analysis, platform)
#            evals.append(os.path.join(eval_path, f'{key}{prefix}results.csv'))
#        print(evals)
#    return pd.concat([pd.read_csv(e) for e in evals if os.path.exists(e)])


def compute_metric(metric, df, include_unmapped=True):

    if not include_unmapped:
        df = df.loc[df.genome != '-', :]

    def precision(df, micro=True):
        cnx_tp = df[df.metric == 'cnx_tp'].value
        cnx_fp = df[df.metric == 'cnx_fp'].value
        if micro:
            total = cnx_tp.sum() + cnx_fp.sum()
            if total == 0:
                return 1
            return cnx_tp.sum() / total
        else:
            total = cnx_tp + cnx_fp
            np.divide(cnx_tp, total, out=np.zeros_like(total, dtype=float), where=total!=0)

    def recall(df, micro=True):
        cnx_tp = df[df.metric == 'cnx_tp'].value
        cnx_fn = df[df.metric == 'cnx_fn'].value
        if micro:
            total = cnx_tp.sum() + cnx_fn.sum()
            if total == 0:
                return 0
            return cnx_tp.sum() / total
        else:
            total = cnx_tp + cnx_fn
            np.divide(cnx_tp, total, out=np.zeros_like(total, dtype=float), where=total!=0)

    def cov_precision(df, micro=True):
        cov_tp = df[df.metric == 'cov_tp'].value
        cov_fp = df[df.metric == 'cov_fp'].value
        if micro:
            total = cov_tp.sum() + cov_fp.sum()
            if total == 0:
                return 1
            return cov_tp.sum() / total
        else:
            total = cov_tp + cov_fp
            np.divide(cov_tp, total, out=np.zeros_like(total, dtype=float), where=total!=0)

    def cov_recall(df, micro=True):
        cov_tp = df[df.metric == 'cov_tp'].value
        cov_fn = df[df.metric == 'cov_fn'].value
        if micro:
            total = cov_tp.sum() + cov_fn.sum()
            if total == 0:
                return 0
            return cov_tp.sum() / total
        else:
            total = cov_tp + cov_fn
            np.divide(cov_tp, total, out=np.zeros_like(total, dtype=float), where=total!=0)

    if metric == 'cnx_F-score':
        p = precision(df)
        r = recall(df)
        return (2*p*r / (p + r)) if (p+r) > 0 else 0
    elif metric == 'cnx_F-score_macro':
        p = precision(df, micro=False)
        r = recall(df, micro=False)
        num = 2*p*r
        den = p+r
        fscore = np.divide(num, den, out=np.zeros_like(den, dtype=float), where=den!=0)
        return np.mean(fscore)

    elif metric == 'cnx_precision':
        return precision(df)
    elif metric == 'cnx_precision_macro':
        return np.mean(precision(df, micro=False))

    elif metric == 'cnx_recall':
        return recall(df)
    elif metric == 'cnx_recall_macro':
        return np.mean(recall(df, micro=False))

    elif metric == 'cov_F-score':
        p = cov_precision(df)
        r = cov_recall(df)
        return (2*p*r / (p + r)) if (p+r) > 0 else 0
    elif metric == 'cov_F-score_macro':
        p = cov_precision(df, micro=False)
        r = cov_recall(df, micro=False)
        num = 2*p*r
        den = p+r
        fscore = np.divide(num, den, out=np.zeros_like(den, dtype=float), where=den!=0)
        return np.mean(fscore)

    elif metric == 'cov_precision':
        return cov_precision(df)
    elif metric == 'cov_precision_macro':
        return np.mean(cov_precision(df, micro=False))

    elif metric == 'cov_recall':
        return cov_recall(df)
    elif metric == 'cov_recall_macro':
        return np.mean(cov_recall(df, micro=False))
    else:
        sys.exit(-1)


def construct_distribution_of_single_metric_across_datasets(metric, raw_data):
    # ensure inputs are valid
    assert metric in [
        'cnx_precision', 'cnx_recall', 'cnx_F-score',
        'cov_precision', 'cov_recall', 'cov_F-score',
    ]

    groupings = ['depth', 'assembler', 'dataset']
    columns = ['depth', 'assembler', 'dataset', 'value']

    # group the table
    groups = raw_data.groupby(groupings)

    # For each group, compute the relevant metric
    scores = pd.DataFrame(index=range(len(groups)), columns=columns)
    for i, (fields, g) in enumerate(groups):
        scores.iloc[i, :] = (*fields, compute_metric(metric, g))

    return scores


def plot_distribution_of_single_metric_across_datasets(data, y_label, fig_name, x_label='assembler', dpi=1400):
    plt.clf()
    ax = sns.boxplot(x=data.assembler, y=data.value, showfliers=False)
    sns.stripplot(x=data.assembler, y=data.value, legend=False, color='black')
    #sns.boxplot(
    #    x=data.assembler, y=data.value,
    #    whiskerprops={'visible': False},
    #    showbox=False, showfliers=False,
    #    showcaps=False, width=0.3,
    #    linewidth=3, ax=ax
    #)


    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.tight_layout()
    plt.savefig(fig_name, dpi=dpi, bbox_inches='tight')
    print(data)
    for a in ['mcvl', 'mh_contigs', 'mh_graph']:
        print(data.loc[data.assembler == 'copangraph', 'value'].values)
        print(data.loc[data.assembler == a, 'value'].values)
        print(wilcoxon(
            x=data.loc[data.assembler == 'copangraph', 'value'].values,
            y=data.loc[data.assembler == a, 'value'].values,
            alternative='greater'
        ))



def graph_quality_plot_depth(dataset, metric, raw_data, group_replicates=False, filter_on=None, hue_val='assembler'):

    # ensure inputs are valid
    plt.clf()
    assert metric in [
        'cnx_precision', 'cnx_recall', 'cnx_F-score',
        'cov_precision', 'cov_recall', 'cov_F-score',
    ]

    if group_replicates: 
        groupings = ['depth', 'assembler']
        columns = ['depth', 'assembler', 'value']
    else:
        groupings = ['depth', 'assembler', 'dataset']
        columns = ['depth', 'assembler', 'dataset', 'value']
    
    # group the table
    groups = raw_data.groupby(groupings)
    
    # For each group, compute the relevant metric
    scores = pd.DataFrame(index=range(len(groups)), columns=columns)
    for i, (fields, g) in enumerate(groups):
        scores.iloc[i, :] = (*fields, compute_metric(metric, g))
    # Plot line graph
    if filter_on:
        scores = scores.loc[scores.assembler == filter_on,:]
    ax = sns.lineplot(x=(scores['depth'].astype(np.float64)), y=scores['value'], hue=scores[hue_val],legend=False)
    ax.set_xlabel('Read pairs')
    ax.set_ylabel(metric)
    ax.set_title(f'{dataset} {metric}')
    #sns.move_legend(ax, 'upper left', bbox_to_anchor=(1, 1))
   
    sum_replicates = '_sum_replicates' if group_replicates else ''
    plt.tight_layout()
    plt.savefig(f'{dataset}_{metric}{sum_replicates}.pdf', dpi=1400, bbox_inches='tight')
    plt.savefig(f'{dataset}_{metric}{sum_replicates}.png', dpi=800, bbox_inches='tight')
    plt.clf()


def graph_quality_plot_coasm(dataset_name, metric, raw_data):

    # ensure inputs are valid
    plt.clf()
    assert metric in [
        'cnx_precision', 'cnx_recall', 'cnx_F-score',
        'cov_precision', 'cov_recall', 'cov_F-score',
    ]

    raw_data.loc[:, 'num_samples'] = raw_data.dataset.apply(lambda x: int(re.findall('([0-9]+)_samples', x)[0]))

    groupings = ['num_samples', 'assembler']
    columns = ['num_samples', 'assembler', 'value']
    
    # group the table
    groups = raw_data.groupby(groupings)
    
    # For each group, compute the relevant metric
    scores = pd.DataFrame(index=range(len(groups)), columns=columns)
    for i, (fields, g) in enumerate(groups):
        scores.iloc[i, :] = (*fields, compute_metric(metric, g))

    # Plot line graph
    ax = sns.lineplot(x=(scores['num_samples'].astype(np.float64)), y=scores['value'], hue=scores['assembler'],legend=True)
    ax.set_xlabel('Number of Samples')
    ax.set_ylabel(metric)
    ax.set_title(f'{dataset_name} {metric}')
    #sns.move_legend(ax, 'upper left', bbox_to_anchor=(1, 1))
   
    plt.tight_layout()
    plt.savefig(f'{dataset_name}_{metric}_coasm.pdf', dpi=1400, bbox_inches='tight')
    plt.savefig(f'{dataset_name}_{metric}_coasm.png', dpi=800, bbox_inches='tight')
    plt.clf()

def build_graph_quality_by_abundance_table(metric, raw_data, group_replicates=False):

    # compute the abundance quantiles. Make quantile bins of 0.1 intervals
    bins = np.linspace(start=0, end=1, num=11)
    quantiles = raw_data.abundance.quantile(q=bins)
    quantiles.index = range(11)

    scores = list()
    # for each assembler...
    for (asm, group) in raw_data.groupby(['assembler']):
        # get the genomes with abundance in a particular quantile, compute the metric, and save in table
        for i in range(10):
            df = group.loc[
                 (quantiles[i] <= raw_data.abundance) & (raw_data.abundance <= quantiles[i+1]), :
            ]
            num_genomes = df.shape[0]
            score_table = pd.DataFrame(
                index=range(num_genomes), columns=['assembler', 'quantile_lower_bound', 'quantile_upper_bound', 'genome', 'value']
            )
            per_genome_metrics = compute_metric(metric, df, per_genome=True)
            score_table.assembler = [asm] * num_genomes
            score_table.quantile_lower_bound = [bins[i]] * num_genomes
            score_table.quantile_upper_bound = [bins[i+1]] * num_genomes
            score_table.value = per_genome_metrics
            scores.append(score_table)

    # concatenate all the scores across assembler and abundance and plot

    return pd.concat(scores)

def plot_graph_quality_by_abundance_table(df, num_assemblies):
    ax = sns.stripplot(x=df.quantile_upper_bound.astype(str), y='value', hue='assembler', data=df, dodge=True, size=2)
    plt.xticks(rotation=45, ha="right")

    # plot the mean line
    sns.boxplot(
        showmeans=True,
        meanline=True,
        meanprops={'color': 'k', 'ls': '-', 'lw': 1.5},
        medianprops={'visible': False},
        whiskerprops={'visible': False},
        zorder=10,
        x="quantile_upper_bound",
        y="value",
        hue='assembler',
        data=df,
        showfliers=False,
        showbox=False,
        showcaps=False,
        ax=ax
    )
    handles, labels = ax.get_legend_handles_labels()
    l = plt.legend(handles[0:num_assemblies], labels[0:num_assemblies]) #, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.xlabel('Relative Abundance')
    plt.ylabel('Metric')
    plt.show()


def concatenate_results_TESTING(results_list):

    with open(results_list) as f:
        data = [pd.read_csv(l.strip()) for l in f]
    return pd.concat(data)


def acu_r01_copan_perfomance_code():
    if len(sys.argv) != 3:
        print('exe <results.csv> <abundance.csv>')
        sys.exit()

    ab = pd.read_csv(sys.argv[2]) if sys.argv[2] != 'none' else None
    above_thresh = ab.loc[ab.abundance > ABD_THRESH, :].header

    csv = pd.read_csv(sys.argv[1])
    csv = csv[csv.genome.isin(list(above_thresh) + ['-'])]

    include = ['COMET-megaHIT-95', 'megaHIT-graph', 'metaCarvel-megaHIT', 'metaSPAdes-graph']
    csv = csv.loc[csv.tool.isin(include), :]
    csv.tool = csv.tool.str.replace('COMET-megaHIT-95', 'COMP graph')
    csv.tool = csv.tool.str.replace('megaHIT-graph', 'megaHIT')
    csv.tool = csv.tool.str.replace('metaSPAdes-graph', 'metaSPAdes')
    csv.tool = csv.tool.str.replace('metaCarvel-megaHIT', 'metaCarvel')
    num_tools = len(csv.tool.value_counts())
    barplot_df = pd.DataFrame(columns=['tool', 'metric', 'value'], index=range(num_tools * 12))

    for i, (tool, g) in enumerate(csv.groupby(['tool'])):
        barplot_df.iloc[12 * i + 0, :] = [tool, 'tp', g[g.metric == 'tp'].value.sum()]
        barplot_df.iloc[12 * i + 1, :] = [tool, 'fp', g[g.metric == 'fp'].value.sum()]
        barplot_df.iloc[12 * i + 2, :] = [tool, 'fn', g[g.metric == 'fn'].value.sum()]
        barplot_df.iloc[12 * i + 3, :] = [tool, 'cov_tp', g[g.metric == 'cov_tp'].value.sum()]
        barplot_df.iloc[12 * i + 4, :] = [tool, 'cov_fp', g[g.metric == 'cov_fp'].value.sum()]
        barplot_df.iloc[12 * i + 5, :] = [tool, 'cov_fn', g[g.metric == 'cov_fn'].value.sum()]
        barplot_df.iloc[12 * i + 6, :] = [
            tool, 'prec', barplot_df.iloc[12 * i + 0, 2] / (barplot_df.iloc[12 * i + 0, 2] + barplot_df.iloc[12 * i + 1, 2])
        ]
        barplot_df.iloc[12 * i + 7, :] = [
            tool, 'recall', barplot_df.iloc[12 * i + 0, 2] / (barplot_df.iloc[12 * i + 0, 2] + barplot_df.iloc[12 * i + 2, 2])
        ]
        p, r = barplot_df.iloc[12 * i + 6, 2], barplot_df.iloc[12 * i + 7, 2]
        barplot_df.iloc[12 * i + 8, :] = [tool, 'fscore', 2*p*r / (p + r)]


        barplot_df.iloc[12 * i + 9, :] = [
            tool, 'cov_prec', barplot_df.iloc[12 * i + 3, 2] / (barplot_df.iloc[12 * i + 3, 2] + barplot_df.iloc[12 * i + 4, 2])
        ]
        barplot_df.iloc[12 * i + 10, :] = [
            tool, 'cov_recall', barplot_df.iloc[12 * i + 3, 2] / (barplot_df.iloc[12 * i + 3, 2] + barplot_df.iloc[12 * i + 5, 2])
        ]
        p, r = barplot_df.iloc[12 * i + 9, 2], barplot_df.iloc[12 * i + 10, 2]
        barplot_df.iloc[12 * i + 11, :] = [tool, 'cov_fscore', 2*p*r / (p + r)]

    #connectivity_filter = (barplot_df.metric == 'fscore')
    #coverage_filter = (barplot_df.metric == 'cov_fscore')

    connectivity_filter = (
            (barplot_df.metric == 'tp') | (barplot_df.metric == 'fp') | (barplot_df.metric == 'fn')
    )
    coverage_filter= (
            (barplot_df.metric == 'cov_tp') | (barplot_df.metric == 'cov_fp') | (barplot_df.metric == 'cov_fn')
    )
    g = sns.catplot(data=barplot_df.loc[connectivity_filter, :], kind='bar', x='tool', y='value', hue='metric')
    g.set_xticklabels(rotation=45, fontsize=8, ha='right')
    ax = g.facet_axis(0,0)
    for container in ax.containers:
        ax.bar_label(container, fontsize=6)
    plt.title(f'abd >= {ABD_THRESH} connectivity')
    g.savefig(f'10G_{ABD_THRESH}_connect.updated.pdf', dpi=1400)
    plt.show()
    plt.clf()

    g = sns.catplot(data=barplot_df.loc[coverage_filter, :], kind='bar', x='tool', y='value', hue='metric')
    g.set_xticklabels(rotation=45, fontsize=8, ha='right')
    ax = g.facet_axis(0,0)
    ax.set(yscale='log')
    #for container in ax.containers:
    #    ax.bar_label(container, fontsize=6, fmt='%.3g')
    plt.title(f'abd >= {ABD_THRESH} coverage')
    g.savefig(f'10G_{ABD_THRESH}_covg.updated.pdf', dpi=1400)
    plt.show()
    plt.clf()


def generate_fake_dataset():
    assemblers = ["asm1", "asm2", "asm3", 'asm4', 'asm5', 'asm6']
    quantiles = [i / 10 for i in range(11)]
    genomes = ["g{}".format(i) for i in range(1, 11)]

    rows = []

    for assembler in assemblers:
        for lower_bound in quantiles:
            for upper_bound in quantiles:
                if lower_bound >= upper_bound:
                    continue

                selected_genomes = random.sample(genomes, 5)
                values = [random.random() for _ in range(5)]

                rows.extend(zip([assembler] * 5, [lower_bound] * 5, [upper_bound] * 5, values, selected_genomes))

    columns = ["assembler", "quantile_lower_bound", "quantile_upper_bound", "value", "genome"]
    df = pd.DataFrame(rows, columns=columns)
    return df


if __name__ == '__main__':

    test_abundance_df = generate_fake_dataset()
    plot_graph_quality_by_abundance_table(test_abundance_df, 6)

