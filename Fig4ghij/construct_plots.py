
import glob
import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

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

ANI = 97
GEN=5000
PREF_COV=f'hmlgy_panel_coverage_ANI{ANI}_gen{GEN}_'
PREF_MULTI=f'hmlgy_panel_multigenome_ANI{ANI}_gen{GEN}_'
OUTDIR='../data/Fig4ghij/outputs'
cov_records = glob.glob('../data/Fig4ghij/*coverage_records.csv')
cov_records = pd.concat(pd.read_csv(e) for e in cov_records)
cov_records.sort_values(by=['rep', 'gen', 'assembler', 'parameter', 'metric'], inplace=True)
cov_records.reset_index(drop=True, inplace=True)
cov_records.gen.value_counts()

subrecord = cov_records.loc[(cov_records.gen == GEN) & (cov_records.assembler == 'megahit') & (cov_records.metric == 'cov_F-score'),:]
subrecord.loc[subrecord.parameter == 141,:]['value'].mean()
subrecord.groupby('parameter')['value'].mean()
plt.figure(figsize=(1.5,2))
ax=sns.lineplot(x=subrecord.parameter, y=subrecord.value, color='red', errorbar='sd')
ax.set_xticks(subrecord.parameter.sort_values())
ax.set_ylim(-.05,1.05)
name = OUTDIR + f'{PREF_COV}_megahit_F-score.pdf'
plt.savefig(name, dpi=1400, bbox_inches='tight')
plot_unlabelled_version(ax, name)
subrecord = cov_records.loc[(cov_records.gen == GEN) & (cov_records.assembler == 'copangraph') & (cov_records.metric == 'cov_F-score'),:]
plt.figure(figsize=(1.5,2))
ax=sns.lineplot(x=subrecord.parameter, y=subrecord.value, color='blue', errorbar='sd')
ax.set_xticks(subrecord.parameter.sort_values())
ax.set_ylim(-.05,1.05)
ax.invert_xaxis()
name = OUTDIR + f'{PREF_COV}_copangraph_F-score.pdf'
plt.savefig(name, dpi=1400, bbox_inches='tight')
plot_unlabelled_version(ax, name)

ms_records = glob.glob('../data/Fig4ghij/*multisample_record.csv')
ms_records = pd.concat(pd.read_csv(e) for e in ms_records)
ms_records.rename({'param':'parameter'}, axis=1, inplace=True)
ms_records.sort_values(by=[ 'assembler', 'parameter'], inplace=True)
ms_records.reset_index(drop=True, inplace=True)
ms_records
ms_records.loc[(ms_records.gen == 9000) & (ms_records.assembler == 'copangraph') & (ms_records.parameter >= 0.05) , :].prop_multigenome_coverage_bp.std()
for GEN in [5000]:
    subrecord = ms_records.loc[(ms_records.gen == GEN)& (ms_records.assembler == 'megahit'), :]
    plt.figure(figsize=(1.5,2))
    ax=sns.lineplot(x=subrecord.parameter, y=subrecord.prop_multigenome_coverage_bp, color='red', errorbar='sd')
    ax.set_xticks(subrecord.parameter.sort_values())
    ax.set_ylim(-.05,1.05)
    name = OUTDIR + f'{PREF_MULTI}_multigenome_coverage_{GEN}_megahit.pdf'
    plt.savefig(name, dpi=1400, bbox_inches='tight')
    plot_unlabelled_version(ax, name)
    
    subrecord = ms_records.loc[(ms_records.gen == GEN)& (ms_records.assembler == 'copangraph'), :]
    plt.figure(figsize=(1.5,2))
    ax=sns.lineplot(x=subrecord.parameter, y=subrecord.prop_multigenome_coverage_bp, color='blue', errorbar='sd')
    ax.set_xticks(subrecord.parameter.sort_values())
    ax.set_ylim(-.05,1.05)
    ax.invert_xaxis()
    name = OUTDIR + f'{PREF_MULTI}_multigenome_coverage_{GEN}_copangraph.pdf'
    plt.savefig(name, dpi=1400, bbox_inches='tight')
    plot_unlabelled_version(ax, name)
