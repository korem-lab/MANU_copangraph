import pandas as pd
import sys
import os
import glob
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


SAMPLES = [30026, 30104, 30526, 30611, 30918, 30939, 30948, 31195, 31280, 31299]
CAMISIMDIR = '../data/Fig3bd/camisim/'
TAXADIR = '../data/Fig4b'

def genome_to_id_tab(s):
    df = pd.read_csv(open(os.path.join(CAMISIMDIR, f'{s}_out/genome_to_id.tsv')), sep='\t', header=None)
    df.columns = ['genome_id', 'file']
    df.file = df.file.apply(lambda x: os.path.basename(x))
    df.index = df.file
    return df

def get_subject_table(fa, g2id):
    file = os.path.basename(fa)
    with open(fa) as f:
        # get all headers
        headers = [l for l in f if l.startswith('>')]
    subjects = [(e.split(' ')[0][1:],file, g2id.loc[file, 'genome_id']) for e in headers]
    return pd.DataFrame(subjects, columns=['subject', 'genome', 'genome_id'])

def parse_taxonomic_profile(s):
    prof = os.path.join(CAMISIMDIR, f'{s}_out/taxonomic_profile_0.txt')
    with open(prof) as f:
        lines = [e.strip() for e in f]
        lines = lines[5:]
        lines = [tuple(e.split('\t')) for e in lines]
    profile = pd.DataFrame(lines, columns=['taxid', 'taxa_rank', 'taxpath', 'taxpathsn', 'percentage', 'genome_id', 'cami_otu'])
    profile = profile.loc[profile.taxa_rank == 'strain', :]
    profile['species'] = profile.taxpathsn.apply(lambda x: x.split('|')[6])
    profile['genera'] = profile.taxpathsn.apply(lambda x: x.split('|')[5])
    profile['family'] = profile.taxpathsn.apply(lambda x: x.split('|')[4])
    profile = profile[['taxid', 'taxa_rank', 'percentage', 'genome_id', 'species', 'genera', 'family']]
    return profile


def get_alns():
    all_alns = list()
    for s in SAMPLES:
        print(s)
        profile = parse_taxonomic_profile(s)
        fastas = glob.glob(os.path.join(CAMISIMDIR, f'{s}_out/genomes/GCA*.fa'))
        g2id = genome_to_id_tab(s)
        fastas = [e for e in fastas if os.path.basename(e) in g2id.index]
        subjects = pd.concat([get_subject_table(e, g2id) for e in fastas])
        subjects = pd.merge(subjects, profile, on='genome_id')
        alns = pd.read_csv(os.path.join(TAXADIR, f'{s}_20M.blast.txt'), sep='\t', header=None)
        alns.columns = ['contig_name', 'subject', 'identity', 'aln_len', 'query_len']
        alns = pd.merge(alns, subjects, on='subject')
        alns['key'] = alns.contig_name.apply(lambda x: f'{s}:{x}')
        alns = alns.loc[((alns.aln_len / alns.query_len) > 0.98) & ((alns.aln_len/alns.query_len) < 1.02) & (alns.identity >= 98), :]
        alns['sample'] = s
        all_alns.append(alns)
    
    all_alns = pd.concat(all_alns).reset_index(drop=True)
    all_alns = all_alns.drop_duplicates()
    return all_alns


def compute_alpha_divs(all_alns):
    # parse copangraph
    sds = ['00', '001', '005', '01', '02', '03', '04', '05', '075', '10', '15', '20', '25', '30']
    alpha_divs= list()
    for sd in sds:
        print(sd)
        with open(os.path.join(TAXADIR, f'graph_sd{sd}.fasta')) as f:
            data = [e.strip().split(':') for e in f if e.startswith('>')]
            data = [(node[1:], sample[:-7], contig) for node, _, sample, _, contig, _, _, _ in data]

        copan_seqs = pd.DataFrame(data, columns=['node', 'sample' ,'contig_name'])
        copan_seqs['key'] = copan_seqs.apply(lambda x: x['sample'] + ':' + x['contig_name'], axis=1)
        copan_alns = pd.merge(copan_seqs[['node', 'key']], all_alns, on='key')
        alphas = copan_alns.groupby('node').apply(lambda x: pd.Series([np.log(x.species.nunique()), np.log(x.genera.nunique()), np.log(x.family.nunique())], index=['alpha_div_sp', 'alpha_div_gn', 'alpha_div_fm']))
        alphas['sd'] = sd
        alpha_divs.append(alphas)
    alpha_divs = pd.concat(alpha_divs)
    return alpha_divs


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

if __name__ == '__main__':

    alpha_divs = pd.read_csv(os.path.join(TAXADIR, 'alpha_divs.csv'), dtype={'sd':str})
    all_alns = pd.read_csv(os.path.join(TAXADIR, 'all_alns.csv'), index_col=0)
    plot_dat = pd.melt(alpha_divs.reset_index(), value_vars=['alpha_div_sp', 'alpha_div_gn', 'alpha_div_fm'], id_vars=['node', 'sd'])
    plot_dat.sd = plot_dat.sd.apply(lambda x: float(f'0.{x}'))

    plt.figure(figsize=(4,4))
    ax = sns.lineplot(x=plot_dat.sd, y=plot_dat.value, hue=plot_dat.variable)
    plt.yscale('log')
    name = os.path.join(TAXADIR, 'taxonomy_curve.pdf')
    ax.set_ylabel('alpha-diversity')
    plt.savefig(name, dpi=1400, bbox_inches='tight')
    plot_unlabelled_version(ax, name)

