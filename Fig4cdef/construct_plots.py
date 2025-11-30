import os
import glob
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

ANI = 99
PREF_COV=f'coverage_magani{ANI}_'
PREF_MULTI=f'multigenome_magani{ANI}_'
OUTDIR = '../data/Fig4cdef'

def plot_unlabelled_version(ax, name, tight_layout=True):
    name = os.path.splitext(name)[0]
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    frame1 = plt.gca()
    frame1.legend().set_visible(False)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(f'{name}_UNLBLD.pdf', dpi=1400, bbox_inches='tight')
    plt.savefig(f'{name}_UNLBLD.png', dpi=900, bbox_inches='tight')
    plt.clf()

cov_records = glob.glob('../data/Fig4cdef/*coverage_records.csv')
cov_records = pd.concat(pd.read_csv(e) for e in cov_records)
cov_records.sort_values(by=[ 'assembler', 'parameter', 'metric'], inplace=True)
cov_records = cov_records.loc[cov_records.parameter != 0,:]
cov_records.reset_index(drop=True, inplace=True)

# F-score by n50
subrecord = cov_records.loc[cov_records.metric == 'cov_F-score', :]
subrecord.to_csv(OUTDIR + PREF_COV + '_coverage_F-score_raw_dat.csv', index=None)
plt.figure(figsize=(4,4))
ax = sns.scatterplot(x=subrecord.n50, y=subrecord.value, hue=subrecord.assembler, s=20)
texts = list()
for i in subrecord.index:
    if subrecord.loc[i, 'assembler'] == 'megahit':
        ano = f'k={int(subrecord.parameter[i])}'
    else:
        ano = f'sd={float(subrecord.parameter[i])}'
    plt.annotate(ano, (subrecord.n50[i], subrecord.value[i]), textcoords='offset points', xytext=(2, 2), size=8)
name = OUTDIR + PREF_COV + 'cov_F-score_n50.pdf'
plt.savefig(name, dpi=1400, bbox_inches='tight')
plot_unlabelled_version(ax, name, tight_layout=False)

# F-score by num_edges
subrecord = cov_records.loc[cov_records.metric == 'cov_F-score', :]
subrecord
plt.figure(figsize=(4,4))
ax = sns.scatterplot(x=subrecord.num_edges, y=subrecord.value, hue=subrecord.assembler, s=20)
ax.set_xscale('log')
texts = list()
for i in subrecord.index:
    if subrecord.loc[i, 'assembler'] == 'megahit':
        ano = f'k={int(subrecord.parameter[i])}'
    else:
        ano = f'sd={float(subrecord.parameter[i])}'
    plt.annotate(ano, (subrecord.num_edges[i], subrecord.value[i]), textcoords='offset points', xytext=(5, 5), size=8)
name = OUTDIR + PREF_COV + 'cov_F-score_num_edges.pdf'
plt.savefig(name, dpi=1400, bbox_inches='tight')
plot_unlabelled_version(ax, name, tight_layout=False)

print('reading mix records')
ms_records = glob.glob('../data/Fig4cdef/*multigenome_record.csv')
ms_records = pd.concat(pd.read_csv(e) for e in ms_records)
ms_records.rename({'param':'parameter'}, axis=1, inplace=True)
ms_records.sort_values(by=[ 'assembler', 'parameter'], inplace=True)
ms_records.reset_index(drop=True, inplace=True)

plt.figure(figsize=(4,4))
ax = sns.scatterplot(x=ms_records.num_edges, y=ms_records.prop_multigenome_coverage_bp, hue=ms_records.assembler, s=20)
ax.set_xscale('log')
texts = list()
for i in ms_records.index:
    if ms_records.loc[i, 'assembler'] == 'megahit':
        ano = f'k={int(ms_records.parameter[i])}'
    else:
        ano = f'sd={float(ms_records.parameter[i])}'
    plt.annotate(ano, (ms_records.num_edges[i], ms_records.prop_multigenome_coverage_bp[i]), textcoords='offset points', xytext=(5, 5), size=8)
name = OUTDIR + PREF_MULTI + 'multi_genome_coverage_edges.pdf'
plt.savefig(name, dpi=1400, bbox_inches='tight')
plot_unlabelled_version(ax, name, tight_layout=False)

plt.figure(figsize=(4,4))
ax = sns.scatterplot(x=ms_records.n50, y=ms_records.prop_multigenome_coverage_bp, hue=ms_records.assembler, s=20)
texts = list()
for i in ms_records.index:
    if ms_records.loc[i, 'assembler'] == 'megahit':
        ano = f'k={int(ms_records.parameter[i])}'
    else:
        ano = f'sd={float(ms_records.parameter[i])}'
    plt.annotate(ano, (ms_records.n50[i], ms_records.prop_multigenome_coverage_bp[i]), textcoords='offset points', xytext=(5, 5), size=8)
name = OUTDIR + PREF_MULTI + 'multi_genome_coverage_n50.pdf'
plt.savefig(name, dpi=1400, bbox_inches='tight')
plot_unlabelled_version(ax, name, tight_layout=False)
