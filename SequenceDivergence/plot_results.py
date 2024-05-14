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


def build_btwn_aln_plot_df(files):

    def add_record(final_df, i, r, sid, asm, depth, dt, is_filtered):

        # mean
        val = np.mean(r)
        final_df.loc[i, :] = [sid, depth, asm, dt, is_filtered, 'mean', pd.NA, val]
        i += 1

        # median
        val = np.median(r)
        final_df.loc[i, :] = [sid, depth, asm, dt, is_filtered, 'median', pd.NA, val]
        i += 1

        # percentile
        for perc in [50, 75, 90, 99, 99.5]:
            val = np.percentile(r, perc)
            final_df.loc[i, :] = [sid, depth, asm, dt, is_filtered, 'percentile', str(perc), val]
            i += 1

        return final_df, i

    final_df = pd.DataFrame(
        index=range(len(files)),
        columns=['sample_id', 'depth', 'asm', 'dt', 'is_filtered', 'summary', 'percentile', 'value']
    )
    idx = 0
    for f in files:
        record = pd.read_pickle(f)
        if len(record) == 10:
            print('skipped', f)
            continue
        sample_id, asm, depth, dt = re.findall('sample_([0-9]+)_(megahit|metaspades)_([0-9]+)M_.*_(0\.[0-9]+).pkl', f)[0]
        sample_id = int(sample_id)
        dt = float(dt)

        # add the non-filtered record
        final_df, idx = add_record(final_df, idx, record, sample_id, asm, depth, dt, False)

        # add the filtered record
        final_df, idx = add_record(final_df, idx, record[record != 1], sample_id, asm, depth, dt, True)

    return final_df


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

        # fill sum, mean, median
        template.summary = 'sum'
        template.non_zero_filtered = False
        template.value = data.sum(axis=1)
        summary_stats.append(template.copy())

        template.summary = 'mean'
        template.non_zero_filtered = False
        template.value = data.mean(axis=1)
        summary_stats.append(template.copy())

        template.summary = 'median'
        template.non_zero_filtered = False
        template.value = data.median(axis=1)
        summary_stats.append(template.copy())

        template.summary = 'mean'
        template.non_zero_filtered = True
        template.value = data.apply(lambda x: x[x!=0].mean(), axis=1)
        summary_stats.append(template.copy())

        template.summary = 'median'
        template.non_zero_filtered = True
        template.value = data.apply(lambda x: x[x!=0].median(), axis=1)
        summary_stats.append(template.copy())

    final_df = pd.concat(summary_stats)
    final_df.index = range(final_df.shape[0])
    return final_df


def plot_aln_percentile_boxplot_for_dt_range_between(files, plot_type, dts, asm):
    plt.clf()
    df = pd.DataFrame(columns=['sample', 'dt', 'percent_above_dt'])
    dts = [str(e) for e in sorted(dts)]
    for dt in sorted(dts):
        for s in range(0, 10):
            file = [e for e in files if dt in e and asm in e and f'sample_{s}' in e]
            if file:
                file = file.pop()
            else:
                continue
            data = pd.read_pickle(file)
            percent_above = 100 - percentileofscore(data, 1-float(dt))
            df.loc[df.shape[0], :] = [s, str(1-float(dt)), percent_above]
    figure(figsize=(4,4), dpi=600)
    sns.boxplot(x=df.dt , y=df.percent_above_dt, showfliers=False)
    sns.stripplot(x=df.dt, y=df.percent_above_dt, jitter=True, color='black')
    plt.xlabel('% threshold')
    plt.ylabel('% sequence above threshold')
    plt.savefig(f'{plot_type}_{asm}_percentile_above_boxplot.pdf')
    plt.tight_layout()
    plt.clf()


def plot_aln_percentile_boxplot_for_dt_range_within(files, plot_type, dts, asm):
    plt.clf()
    df = pd.DataFrame(columns=['sample', 'dt', 'percent_above_dt'])
    dts = [str(e) for e in sorted(dts)]
    for dt in sorted(dts):
        for s in range(0, 10):
            file = [e for e in files if dt in e and asm in e and f'sample_{s}' in e]
            if file:
                file = file.pop()
            else:
                continue
            r = pd.read_pickle(file).iloc[0,:]
            filt = [type(e) == int for e in r.index]
            data = r[filt]
            data = data[(data != 0) & (data != np.inf)]
            percent_above = 100 - percentileofscore(data, 1-float(dt))
            df.loc[df.shape[0], :] = [s, str(1-float(dt)), percent_above]
    figure(figsize=(3,3), dpi=600)
    plt.xlabel('% threshold')
    plt.ylabel('% sequence above threshold')
    sns.boxplot(x=df.dt , y=df.percent_above_dt, showfliers=False)
    sns.stripplot(x=df.dt, y=df.percent_above_dt, jitter=True, color='black')
    plt.savefig(f'{plot_type}_{asm}_percentile_above_boxplot.pdf')
    plt.tight_layout()
    plt.clf()


def plot_aln_distribution_between(files, plot_type, dt, asm, min_x=0.5):
    files = [e for e in files if str(dt) in e and asm in e]
    print(files)
    dt = 1-dt
    plt.clf()
    figure(figsize=(6,2), dpi=600)
    for f in files:
        data = pd.read_pickle(f)
        kde = sns.kdeplot(data)
        kde_x, kde_y = kde.get_lines()[0].get_data()
        plt.fill_between(kde_x, kde_y, where=(kde_x >= dt), interpolate = True, color='darkblue', alpha=0.50)
        plt.axvline(dt, color='darkred', linestyle='--')
        plt.xlim((min_x, 1.0))
        plt.xlabel(f'{dt}% threshold')
        plt.tight_layout()
        sample_id = re.findall('sample_([0-9])_', f)[0]
        plt.savefig(f'sample_{sample_id}_{plot_type}_{dt}_{asm}.pdf')
        plt.clf()


def plot_aln_distribution_within(files, plot_type, dt, asm, min_x=0.95):
    files = [e for e in files if str(dt) in e and asm in e]
    records = [pd.read_pickle(e).iloc[0, :] for e in files]
    dt = 1 - dt
    plt.clf()
    figure(figsize=(6, 2), dpi=600)
    for r in records:
        # replace sample id with numeric value
        r.sample_id = int(re.findall('sample_([0-9]+)_', r.sample_id)[0])
        filt = [type(e) == int for e in r.index]
        data = r[filt]
        data = data[(data != 0) & (data != np.inf) & (data >= min_x)]
        kde = sns.kdeplot(data)
        kde_x, kde_y = kde.get_lines()[0].get_data()
        plt.fill_between(kde_x, kde_y, where=(kde_x >= dt), interpolate = True, color='darkblue', alpha=0.50)
        plt.axvline(dt, color='darkred', linestyle='--')
        plt.xlim((min_x, data.max()))
        plt.xlabel(f'{dt}% threshold')
        plt.tight_layout()
        plt.savefig(f'sample_{r.sample_id}_{plot_type}_{dt}_{asm}.pdf')
        plt.clf()


def build_aln_plot_df(files):
    # all columns are identical, so pick the first
    records = [pd.read_pickle(e).iloc[0, :] for e in files]
    summary_stats = list()
    for r in records:

        # replace sample id with numeric value
        r.sample_id = int(re.findall('sample_([0-9]+)_', r.sample_id)[0])

        # construct summary stats
        template = pd.Series(index=['sample_id', 'depth', 'asm', 'dt', 'summary', 'percentileofscore', 'percentile', 'value'], dtype=object)
        template.sample_id = r.sample_id
        template.depth = r.depth
        template.asm = r.asm
        template['dt'] = r['dt']
        # extract just the alignment scores
        filt = [type(e) == int for e in r.index]
        data = r[filt]
        data = data[(data != 0) & (data != np.inf)]   # filter zero or inf values which represent nodes with a single sequence

        # construct mean, filter zero values which represent nodes with a single sequence
        template.value = data.mean()
        template.summary = 'mean'
        summary_stats.append(template.copy())

        # construct median
        template.value = data.median()
        template.summary = 'median'
        summary_stats.append(template.copy())

        # construct percentileofscore
        for pos in [0.0, 0.005, 0.01, 0.02]:
            template.value = (100 - percentileofscore(data, 1-(template['dt'] + pos)))/100.0
            template.percentileofscore = f'dt + {pos}'
            template.summary = 'percentileofscore'
            summary_stats.append(template.copy())


        # construct percentiles
        for perc in [8, 10, 25, 50]:
            template.value = np.percentile(data, perc)
            template.summary = 'percentile'
            template.percentile = str(perc)
            summary_stats.append(template.copy())

    final_df = pd.DataFrame(summary_stats)
    return final_df


def plot_btwn_figs(btwn_plot_df, asm='megahit', is_filtered=False):
    assert asm in ['megahit', 'metaspades']
    assert is_filtered in [True, False]

    filtered = 'exclude_exact' if is_filtered else 'include_exact'

    # filter for assembler
    btwn_plot_df = btwn_plot_df.loc[(btwn_plot_df.asm == asm) & (btwn_plot_df.is_filtered == is_filtered), :]

    # plot mean
    dat = btwn_plot_df.loc[btwn_plot_df.summary == 'mean', :]
    ax = sns.lineplot(x=1-dat['dt'], y=dat.value)
    plt.ylim(0, 1)
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(f'btwn_seqid_mean_{asm}_{filtered}.pdf')
    plt.clf()

    # plot median
    dat = btwn_plot_df.loc[btwn_plot_df.summary == 'median', :]
    ax = sns.lineplot(x=1-dat['dt'], y=dat.value)
    plt.ylim(0, 1)
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(f'btwn_seqid_median_{asm}_{filtered}.pdf')
    plt.clf()

    # plot percentile
    dat = btwn_plot_df.loc[btwn_plot_df.summary == 'percentile', :]
    ax = sns.lineplot(x=1-dat['dt'], y=dat.value, hue=dat.percentile)
    plt.ylim(0, 1)
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(f'btwn_seqid_percentile_{asm}_{filtered}.pdf')
    plt.clf()


def plot_figs(plot_df, plot_type, asm='megahit', non_zero_filtered=False):

    assert plot_type in ['entropies', 'ani_within', 'ani_btwn']
    assert asm in ['megahit', 'metaspades']

    # filter for assembler
    plot_df = plot_df.loc[plot_df.asm == asm, :]
    plot_df = plot_df.loc[(plot_df.tax_lvl == 'Species') | (plot_df.tax_lvl == 'Genera') | (plot_df.tax_lvl == 'Family'), :]

    # filter for zero
    plot_df = plot_df.loc[plot_df.non_zero_filtered == non_zero_filtered, :]

    filtered = 'non_zero_filtered' if non_zero_filtered else 'unfiltered'
    plt.clf()
    figure(figsize=(4,4),dpi=600)
    # plot sum
    if not non_zero_filtered:
        dat = plot_df.loc[plot_df.summary == 'sum', :]
        ax = sns.lineplot(x=dat['dt'], y=dat.value, hue=dat.tax_lvl)
        plt.tight_layout()
        plt.savefig(f'{plot_type}_sum_{asm}_{filtered}.pdf')
        plt.clf()
    # plot mean
    dat = plot_df.loc[plot_df.summary == 'mean', :]
    ax = sns.lineplot(x=dat['dt'], y=dat.value, hue=dat.tax_lvl,
                      estimator='median', errorbar=lambda x: (np.percentile(x, 10), np.percentile(x, 90)), legend=True)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(f'{plot_type}_mean_{asm}_{filtered}_leg_plotlog10.pdf')
    plt.clf()

    # plot median
    dat = plot_df.loc[plot_df.summary == 'median', :]
    ax = sns.lineplot(x=dat['dt'], y=dat.value, hue=dat.tax_lvl)
    plt.tight_layout()
    plt.savefig(f'{plot_type}_median_{asm}_{filtered}.pdf')
    plt.clf()


def entropies_boxplot(files, dts, tax_levels, asm):

    df = pd.DataFrame(columns=['sample', 'dt', 'tax_level', 'mean_entropy'])
    files = sorted(files)
    for dt in dts:
        dt_files = [e for e in files if f'{dt}.pkl' in e and asm in e]
        data = [pd.read_pickle(e) for e in dt_files]
        for r in data:
            sample_id = int(re.findall('sample_([0-9]+)', r.sample_id[0])[0])
            filt = [type(e) == int for e in r.columns]
            data = r.loc[:, filt]
            for tax in tax_levels:
                df.loc[df.shape[0], :] = [sample_id, str(1 - r.dt[0]), tax, data.loc[tax, :].mean()]
    plt.clf()
    figure(figsize=(6,3), dpi=600)
    dt_pairs = zip(dts, dts[1:])
    for tax in tax_levels:
        df_tax = df.loc[df.tax_level == tax,:]
        sns.boxplot(x=df_tax.dt, y=df_tax.mean_entropy)
        sns.stripplot(x=df_tax.dt, y=df_tax.mean_entropy, jitter=True, color='black')
        plt.tight_layout()
        plt.savefig(f'entropies_boxplot_{tax}_{asm}.pdf')
        plt.clf()
        for a, b in dt_pairs:
            df_tax_x = df_tax.loc[(df_tax.dt.astype(str) == str(1-a)), ]
            df_tax_y = df_tax.loc[(df_tax.dt.astype(str) == str(1-b)), ]
            df_tax_merge = pd.merge(df_tax_x, df_tax_y, on='sample')
            print(df_tax_x, df_tax_y, df_tax_merge)
            sys.exit()
            #print(
            #    tax, str(1-a), str(1-b), scipy.stats.wilcoxon(
            #        x=df_tax.loc[df_tax.dt.astype(str) == str(1-a), 'mean_entropy'],
            #        y=df_tax.loc[df_tax.dt.astype(str) == str(1-b), 'mean_entropy']
            #    )
            #)




def plot_aln_figs(aln_plot_df, plot_type, asm='megahit'):

    assert plot_type in ['aln_avpv', 'aln_minpv']
    assert asm in ['megahit', 'metaspades']

    # filter for assembler
    aln_plot_df = aln_plot_df.loc[aln_plot_df.asm == asm, :]

    # plot mean
    dat = aln_plot_df.loc[aln_plot_df.summary == 'mean', :]
    ax = sns.lineplot(x= 1- dat['dt'], y=dat.value)
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(f'{plot_type}_mean_{asm}.pdf')
    plt.clf()

    # plot median
    dat = aln_plot_df.loc[aln_plot_df.summary == 'median', :]
    ax = sns.lineplot(x=1- dat['dt'], y=dat.value)
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(f'{plot_type}_median_{asm}.pdf')
    plt.clf()

    # plot percentileofscore
    dat = aln_plot_df.loc[aln_plot_df.summary == 'percentileofscore', :]
    ax = sns.lineplot(x=1-dat['dt'], y=dat.value, hue=dat.percentileofscore)
    plt.ylim(0, 1)
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(f'{plot_type}_percentileofscore_{asm}.pdf')
    plt.clf()

    # plot percentile
    dat = aln_plot_df.loc[aln_plot_df.summary == 'percentile', :]
    ax = sns.lineplot(x=1-dat['dt'], y=dat.value, hue=dat.percentile)
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(f'{plot_type}_percentile_{asm}.pdf')
    plt.clf()


if __name__ == '__main__':

    if len(sys.argv) != 3:
        print('exe: <out_dir> <depth>')
        sys.exit()

    _, out_dir, depth = sys.argv
    all_files = glob.glob(os.path.join(out_dir, '*.pkl'))
    all_files = [e for e in all_files if depth in e]

    ## plot aln_avpv
    print('plotting aln_avpv')
    plot_type = 'aln_avpv'
    ##aln_plot_df = build_aln_plot_df([e for e in all_files if plot_type in e])
    ##aln_plot_df.to_csv(f'{plot_type}_plot_df.csv', index=None)
    ##plot_aln_figs(aln_plot_df, plot_type, asm='megahit')
    ##plot_aln_figs(aln_plot_df, plot_type, asm='metaspades')
    #plot_aln_percentile_boxplot_for_dt_range_within([e for e in all_files if plot_type in e], plot_type, [0.01, 0.03, 0.05], 'megahit')
    #plot_aln_percentile_boxplot_for_dt_range_within([e for e in all_files if plot_type in e], plot_type, [0.01, 0.03, 0.05], 'metaspades')
    #plot_aln_distribution([e for e in all_files if plot_type in e], plot_type, dt=0.01, asm='megahit', min_x=0.9)
    #plot_aln_distribution([e for e in all_files if plot_type in e], plot_type, dt=0.03, asm='megahit', min_x=0.9)
    #plot_aln_distribution([e for e in all_files if plot_type in e], plot_type, dt=0.05, asm='megahit', min_x=0.9)
    #plot_aln_distribution([e for e in all_files if plot_type in e], plot_type, dt=0.01, asm='metaspades', min_x=0.9)
    #plot_aln_distribution([e for e in all_files if plot_type in e], plot_type, dt=0.03, asm='metaspades', min_x=0.9)
    #plot_aln_distribution([e for e in all_files if plot_type in e], plot_type, dt=0.05, asm='metaspades', min_x=0.9)

    #### plot aln_minpv
    print('plotting aln_minpv')
    plot_type = 'aln_minpv'
    ##aln_plot_df = build_aln_plot_df([e for e in all_files if plot_type in e])
    ##aln_plot_df.to_csv(f'{plot_type}_plot_df.csv', index=None)
    ##plot_aln_figs(aln_plot_df, plot_type, asm='megahit')
    ##plot_aln_figs(aln_plot_df, plot_type, asm='metaspades')
    #plot_aln_percentile_boxplot_for_dt_range_within([e for e in all_files if plot_type in e], plot_type, [0.01, 0.03, 0.05], 'megahit')
    #plot_aln_percentile_boxplot_for_dt_range_within([e for e in all_files if plot_type in e], plot_type, [0.01, 0.03, 0.05], 'metaspades')
    #plot_aln_distribution_within([e for e in all_files if plot_type in e], plot_type, dt=0.01, asm='megahit', min_x=0.9)
    #plot_aln_distribution_within([e for e in all_files if plot_type in e], plot_type, dt=0.03, asm='megahit', min_x=0.9)
    #plot_aln_distribution_within([e for e in all_files if plot_type in e], plot_type, dt=0.05, asm='megahit', min_x=0.9)
    #plot_aln_distribution_within([e for e in all_files if plot_type in e], plot_type, dt=0.01, asm='metaspades', min_x=0.9)
    #plot_aln_distribution_within([e for e in all_files if plot_type in e], plot_type, dt=0.03, asm='metaspades', min_x=0.9)
    #plot_aln_distribution_within([e for e in all_files if plot_type in e], plot_type, dt=0.05, asm='metaspades', min_x=0.9)

    ## plot btwn_seqid
    print('plotting btwn_seqid')
    plot_type = 'btwn_seqid'
    #btwn_plot_df = build_btwn_aln_plot_df([e for e in all_files if plot_type in e])
    #btwn_plot_df.to_csv(f'{plot_type}_plot_df.csv')
    #plot_aln_percentile_boxplot_for_dt_range_between([e for e in all_files if plot_type in e], plot_type, [0.01, 0.03, 0.05], 'megahit')
    #plot_aln_percentile_boxplot_for_dt_range_between([e for e in all_files if plot_type in e], plot_type, [0.01, 0.03, 0.05], 'metaspades')
    #plot_aln_distribution_between([e for e in all_files if plot_type in e], plot_type, dt=0.01, asm='megahit', min_x=0.4)
    #plot_aln_distribution_between([e for e in all_files if plot_type in e], plot_type, dt=0.03, asm='megahit', min_x=0.4)
    #plot_aln_distribution_between([e for e in all_files if plot_type in e], plot_type, dt=0.05, asm='megahit', min_x=0.4)
    #plot_aln_distribution_between([e for e in all_files if plot_type in e], plot_type, dt=0.01, asm='metaspades', min_x=0.4)
    #plot_aln_distribution_between([e for e in all_files if plot_type in e], plot_type, dt=0.03, asm='metaspades', min_x=0.4)
    #plot_aln_distribution_between([e for e in all_files if plot_type in e], plot_type, dt=0.05, asm='metaspades', min_x=0.4)
    #plot_btwn_figs(btwn_plot_df, asm='megahit', is_filtered=False)
    #plot_btwn_figs(btwn_plot_df, asm='metaspades', is_filtered=False)
    #plot_btwn_figs(btwn_plot_df, asm='megahit', is_filtered=True)
    #plot_btwn_figs(btwn_plot_df, asm='metaspades', is_filtered=True)

    ## plot entropies
    #print('plotting entropies')
    plot_type = 'entropies'
    #entropies_boxplot([e for e in all_files if plot_type in e], [0.01, 0.03, 0.05, 0.1], ['Strain', 'Species', 'Genera', 'Family'], asm='megahit')
    plot_df = build_plot_df([e for e in all_files if plot_type in e])
    #btwn_plot_df.to_csv(f'{plot_type}_plot_df.csv')
    plot_figs(plot_df, plot_type, asm='megahit', non_zero_filtered=False)
    plot_figs(plot_df, plot_type, asm='metaspades', non_zero_filtered=False)
    #plot_figs(btwn_plot_df, plot_type, asm='metaspades', non_zero_filtered=False)
    #plot_figs(btwn_plot_df, plot_type, asm='megahit', non_zero_filtered=True)
    #plot_figs(btwn_plot_df, plot_type, asm='metaspades', non_zero_filtered=True)


    ## plot ANI  within
    #print('plotting ani within')
    #plot_type = '10M_ani_avpv'
    #btwn_plot_df = build_plot_df([e for e in all_files if plot_type in e])
    #btwn_plot_df.to_csv(f'ani_within_plot_df.csv')
    #plot_figs(btwn_plot_df, 'ani_within', asm='megahit', non_zero_filtered=False)
    #plot_figs(btwn_plot_df, 'ani_within', asm='metaspades', non_zero_filtered=False)
    #plot_figs(btwn_plot_df, 'ani_within', asm='megahit', non_zero_filtered=True)
    #plot_figs(btwn_plot_df, 'ani_within', asm='metaspades', non_zero_filtered=True)

    ## plot ANI  between
    #print('plotting ani between')
    #plot_type = '10M_btwn_ani_avpv'
    #btwn_plot_df = build_plot_df([e for e in all_files if plot_type in e])
    #btwn_plot_df.to_csv(f'ani_btwn_plot_df.csv')
    #plot_figs(btwn_plot_df, 'ani_btwn', asm='megahit', non_zero_filtered=False)
    #plot_figs(btwn_plot_df, 'ani_btwn', asm='metaspades', non_zero_filtered=False)
    #plot_figs(btwn_plot_df, 'ani_btwn', asm='megahit', non_zero_filtered=True)
    #plot_figs(btwn_plot_df, 'ani_btwn', asm='metaspades', non_zero_filtered=True)


#def plot_data(all_files, out, plot_type, filter_zeros=False, average='mean'):
#    files = [e for e in all_files if plot_type in e]
#    dfs = [pd.read_pickle(e) for e in files]
#    for i, f in enumerate(files):
#        dfs[i].sample_id = dfs[i]['sample_id'].apply(lambda x: re.findall('sample_([0-9]+)', x))
#
#    all_avs = list()
#    for df in dfs:
#        data = df.loc[:, [type(e) == int for e in df.columns]]
#        if average == 'mean':
#            avg = data.apply(lambda x: (x[(x != 0) & (x != np.inf)].mean()) if filter_zeros else x.mean(), axis=1)
#        elif average == 'median':
#            avg = data.apply(lambda x: (x[(x != 0) & (x != np.inf)].median()) if filter_zeros else x.median(), axis=1)
#        else:
#            Exception(f"{average} is not a valid average")
#        avg = pd.DataFrame(avg, columns=[plot_type])
#        avg['dt'] = df['dt']
#        avg['sample_id'] = df['sample_id']
#        all_avs.append(avg)
#    avgs = pd.concat(all_avs)
#    avgs['taxa_lvl'] = avgs.index
#    avgs.index = range(avgs.shape[0])
#    ax = sns.lineplot(x=avgs.dt, y=avgs[plot_type], hue=avgs.taxa_lvl)
#    plt.tight_layout()
#    plt.savefig(os.path.join(out, f'{plot_type}_plot_{average}.pdf'))
#    plt.clf()
#
#def plot_aln(all_files, out, plot_type, average='mean'):
#    files = [e for e in all_files if plot_type in e]
#    dfs = [pd.read_pickle(e) for e in files]
#    for i, f in enumerate(files):
#        dfs[i].sample_id = dfs[i]['sample_id'].apply(lambda x: re.findall('sample_([0-9]+)', x))
#        df[i] = df.loc[0, :]  # all columns are identical in each df, so pick the first
#    all_avs = list()
#    if summary == 'mean' or summary == 'median':
#        avg = pd.DataFrame(index=range(len(dfs)), columns=['plot_type', 'summary', 'dt', 'sample_id'])
#        avg.summary = summary
#        avg.dt = [df.dt for df in dfs]
#        avg.sample_id = [df.sample_id for df in dfs]
#        values = list()
#        for df in dfs:
#            df = df.loc[:, type(e) == int for e in df.columns]
#            val = df[df!=0].mean() if summary == 'mean' else df[df!=0].median()
#            values.append(val)
#        avg.plot_type = values
#
#    for df in dfs:
#        data = df.loc[:, [type(e) == int for e in df.columns]]
#        if summary == 'mean':
#            avg = data.apply(lambda x: (x[(x != 0)].mean()), axis=1)
#        elif summary == 'median':
#            avg = data.apply(lambda x: (x[(x != 0)].median()), axis=1)
#        elif summary == 'percentiles':
#            percs = list()
#            for perc in [5, 25, 50]:
#                percs.append(
#                    data.apply(lambda x: np.perc(x[(x != 0)], perc), axis=1)
#                )
#            avg = pd.concat(percs)
#
#
#        avg = pd.DataFrame(avg, columns=[plot_type])
#        avg['dt'] = df['dt']
#        avg['sample_id'] = df['sample_id']
#        all_avs.append(avg)
#    avgs = pd.concat(all_avs)
#    avgs.index = range(avgs.shape[0])
#    ax = sns.lineplot(x=avgs.dt, y=avgs[plot_type])
#    plt.tight_layout()
#    plt.savefig(os.path.join(out, f'{plot_type}_plot_{average}.pdf'))
#    plt.clf()
#
#def distribution(all_files, out, plot_type, filter_zeros=False):
#    files = [e for e in all_files if plot_type in e]
#    dfs = [pd.read_pickle(e) for e in files]
#    for i, f in enumerate(files):
#        dt = re.findall('_(0\.[0-9]+)\.pkl', f)[0]
#        dt = float(dt)
#        dfs[i]['dt'] = dt
#    dfs = sorted(dfs, key=lambda x: x['dt'][0])
#    data = [df.loc[:, [type(e) == int for e in df.columns]] for df in dfs]
#    data = [df.loc[df.index == 'Kingdom', :] for df in data]
#    data = [df.loc[:, df.iloc[0,:] != 0] for df in data]
#    percentiles = [(0, 100), (1, 99), (10, 90), (25, 50)]
#    median = 50
#    divergence = np.zeros(len(data))
#    for j, (l, u) in enumerate(reversed(percentiles), start=1):
#        lower_percentiles = np.zeros(len(data))
#        upper_percentiles = np.zeros(len(data))
#        for i, d in enumerate(data):
#            lower_percentiles[i] = np.percentile(d, q=l)
#            upper_percentiles[i] = np.percentile(d, q=u)
#            divergence[i] = dfs[i]['dt'][0]
#
#        plt.fill_between(
#            divergence, upper_percentiles, lower_percentiles,
#            where=(upper_percentiles > lower_percentiles), interpolate=True,
#            color='green', alpha=(float(j) * 0.3/ (len(percentiles)+1)), edgecolor='none'
#        )
#    print(divergence)
#    median_percentile=np.zeros(len(data))
#    for i, d in enumerate(data):
#        median_percentile[i] = np.percentile(d, q=median)
#    plt.plot(divergence, median_percentile, color='green')
#    plt.savefig('testing.pdf')
#
#
#
#if __name__ == '__main__':
#
#    if len(sys.argv) != 4:
#        print('exe: <plot_type> <results_dir> <out>')
#        sys.exit()
#
#    _, plot_type, result_dir, out_dir = sys.argv
#
#    all_files = glob.glob(os.path.join(result_dir, '*.pkl'))
#    if plot_type == 'entropies':
#        plot_data(all_files, out_dir, plot_type, average='mean')
#        plot_data(all_files, out_dir, plot_type, average='median')
#    elif plot_type == 'ani_avpv':
#        plot_data(all_files, out_dir, plot_type, average='mean')
#        plot_data(all_files, out_dir, plot_type, average='median')
#    elif plot_type == 'aln_avpv':
#        plot_data(all_files, out_dir, plot_type, filter_zeros=True, average='mean')
#        plot_data(all_files, out_dir, plot_type, filter_zeros=True, average='median')
#        #distribution(all_files, out_dir, plot_type)
#    elif plot_type == 'aln_minpv':
#        plot_data(all_files, out_dir, plot_type, filter_zeros=True, average='mean')
#        plot_data(all_files, out_dir, plot_type, filter_zeros=True, average='median')
#        #distribution(all_files, out_dir, plot_type)




