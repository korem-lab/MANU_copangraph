import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import tqdm
import glob
from scipy.spatial import distance as dist
import os
from uuid import uuid4
from subprocess import run
from skbio.stats.ordination import pcoa
from matplotlib.ticker import MaxNLocator



def read_copan(file, index):
    if '.pkl' in file:
        df = pd.read_pickle(file).T
    else:
        cols = pd.read_csv(file, nrows=0).columns[1:]
        dtype_map = {col: np.int8 for col in cols}
        df = pd.read_csv(file, index_col=0, dtype=dtype_map).T
    return df.loc[index,:]

def remove_correlated_at(X, corr=1):
    assert corr == 1
    if corr == 1:
        X = X.loc[:, ~X.T.duplicated()]
    return X

def occ_filter(df, low_perc, high_perc, t, pref):
    if 'abnd' not in pref:
        t=0
    print(pref, t)
    low_perc /=100.0
    high_perc /=100.0
    occ = (df>t).sum()
    df = df.loc[:, (occ >= df.shape[0]*low_perc) & (occ <= df.shape[0]*high_perc)]
    return df

def is_constant(x):
    return np.allclose(x, x[0])

def get_dist(df, distance):
    if distance == 'euclidean':
        return dist.cdist(df.values, df.values, distance)
    elif distance == 'braycurtis':
        D = dist.cdist(df.values + 1e-6, df.values + 1e-6, distance)
        return np.clip((1-np.eye(D.shape[0]))*D, 0,1)
    else:
        D = dist.cdist(df.values, df.values, distance)
        return np.clip((1-np.eye(D.shape[0]))*D, 0,1)
    
def run_adonis(X_fl, outcome_fl, is_factor, is_dist, X_nm = '../data/Fig6abcdef/adonis/X_', y_nm='../data/Fig6abcdef/adonis/y_'):
    os.makedirs('../data/Fig6abcdef/adonis',exist_ok=True)
    _uuid=str(uuid4())
    X_nm+=_uuid +'.csv'
    X_fl.to_csv(X_nm)
    y_nm+=_uuid +'.csv'
    outcome_fl.to_csv(y_nm)
    cmd = f'Rscript adonis.R {X_nm} {y_nm} {"factor" if is_factor else "not_factor"} {"dist" if is_dist else "not_dist" }'#> /dev/null 2>&1'
    print(cmd)
    run(cmd, shell=True)
    adonis_fl = X_nm.replace('.csv', '_adonis.csv')
    os.remove(X_nm)
    os.remove(y_nm)
    adonis = pd.read_csv(adonis_fl)
    os.remove(adonis_fl)
    return adonis

def plot_unlabelled_version(ax, name, tight_layout=True):
    name = os.path.splitext(name)[0]
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')
    frame1 = plt.gca()
    frame1.legend().set_visible(False)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(f'{name}_UNLBLD.pdf', dpi=1400, bbox_inches='tight')
    plt.savefig(f'{name}_UNLBLD.png', dpi=900, bbox_inches='tight')
    plt.clf()


if __name__ == '__main__':

    #EUCLID DISTANCE TOP CONFIG
    labels = pd.read_csv('../data/Fig6abcdef/momspi/Labels.csv', index_col=0)
    df = read_copan('../data/Fig6abcdef/momspi/graph_euclid.abnd_mat.pkl', labels.index.astype(str))
    df=occ_filter(df,20,33,50,'.abnd_mat.pkl')
    D = get_dist(df,'euclidean')
    adonis =run_adonis(pd.DataFrame(D, index=labels.index, columns=labels.index), labels, True, True) 
    ordination = pcoa(D)
    coords = ordination.samples
    coords['outcome'] = labels.outcome.values
    plt.figure(figsize=(2.5,2.5))
    ax = sns.scatterplot(data=coords, x='PC1', y='PC2', hue='outcome',alpha=0.4,s=50)
    ax.set_title(f'mpi e r:{adonis.R2[0]:.4f},F:{adonis.F[0]:.4f},p:{adonis["Pr(>F)"][0]},N:{df.shape[1]}')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    name='../data/Fig6abcdef/Fig6f'
    print('PC1 - ', ordination.proportion_explained['PC1'], 'PC2 - ', ordination.proportion_explained['PC2'])
    plt.savefig(name+'_LBL.pdf',dpi=1400, bbox_inches='tight')
    plt.savefig(name+'_LBL.png',dpi=900, bbox_inches='tight')
    plot_unlabelled_version(ax,name)

    #JACCARD DISTANCE TOP CONFIG
    df = read_copan('../data/Fig6abcdef/momspi/graph_jaccard.abnd_mat.pkl', labels.index.astype(str))
    df=occ_filter(df,20,33,50,'.abnd_mat.pkl')
    D = get_dist(df,'jaccard')
    adonis =run_adonis(pd.DataFrame(D, index=labels.index, columns=labels.index), labels, True, True) 
    ordination = pcoa(D)
    coords = ordination.samples
    coords['outcome'] = labels.outcome.values
    plt.figure(figsize=(2.5,2.5))
    ax = sns.scatterplot(data=coords, x='PC1', y='PC2', hue='outcome',alpha=0.4,s=50)
    ax.set_title(f'mpi j r:{adonis.R2[0]:.4f},F:{adonis.F[0]:.4f},p:{adonis["Pr(>F)"][0]},N:{df.shape[1]}')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    name='../data/Fig6abcdef/Fig6e'
    print('PC1 - ', ordination.proportion_explained['PC1'], 'PC2 - ', ordination.proportion_explained['PC2'])
    plt.savefig(name+'_LBL.pdf',dpi=1400, bbox_inches='tight')
    plt.savefig(name+'_LBL.png',dpi=900, bbox_inches='tight')
    plot_unlabelled_version(ax,name)

    df = pd.read_pickle('../data/Fig6abcdef/momspi/mOTUs.df')
    df=df.loc[labels.index,:]
    D=get_dist(df,'braycurtis')
    adonis=run_adonis(pd.DataFrame(D, index=labels.index, columns=labels.index), labels, True, True) 
    ordination = pcoa(D)
    coords= ordination.samples
    coords['outcome'] = labels.outcome.values
    print('PC1 - ', ordination.proportion_explained['PC1'], 'PC2 - ', ordination.proportion_explained['PC2'])
    plt.figure(figsize=(2.5,2.5))
    ax = sns.scatterplot(data=coords, x='PC1', y='PC2', hue='outcome',alpha=0.4,s=50)
    ax.set_title(f'mpi bc r:{adonis.R2[0]:.4f},F:{adonis.F[0]:.4f},p:{adonis["Pr(>F)"][0]},N:{df.shape[1]}')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    name='../data/Fig6abcdef/Fig6d'
    plt.savefig(name+'_LBL.pdf',dpi=1400, bbox_inches='tight')
    plt.savefig(name+'_LBL.png',dpi=900, bbox_inches='tight')
    plot_unlabelled_version(ax,name)
