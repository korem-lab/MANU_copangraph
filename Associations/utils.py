import scipy.stats as stats
import numpy as np
import seaborn as sns
import pandas as pd
import parse_seq as ps
import tqdm
from statsmodels.stats.multitest import multipletests
from scipy.sparse import lil_matrix
from scipy.cluster.hierarchy import linkage, fcluster

def get_outcome(path):
    if 'moms-pi' in path:
        outcome = pd.read_csv(path)
        outcome.index = outcome['sample']
        outcome = outcome.outcome
        return outcome
    if 'ACU' in path:
        outcome = pd.read_csv(path, index_col=0)
        outcome.index = outcome.index.astype(str)
        outcome = outcome.groupby('StudyID').head(1)  # get one sample per patient
        outcome = outcome['Outcome'].apply(lambda x: 1 if x == 'Persistence' else 0)
        return outcome
    if 'CRC' in path:
        df = pd.read_excel(path)
        outcome = df.config.apply(lambda x: 1 if x == 'case' else 0)
        outcome.index = df.Run
        return outcome
    

def top_percentile_filter(X, perc=99.5):
    """ Excludes columns whose sum is >= perc of all column sums"""
    p = np.percentile(np.sort(X.sum(0)), perc)
    X = X.loc[:, X.sum(0) <= p]
    return X


def get_cluster_averages(X, corr=0.9):
    link_mat = linkage(X.T, method='complete', metric='correlation')
    clusters = fcluster(link_mat, t=1-corr, criterion='distance')
    clust_av = list()
    for c in np.unique(clusters):
        c_df = X.loc[:, clusters == c].mean(1)
        clust_av.append(c_df)
    return pd.concat(clust_av, axis=1)



def abundance_filter(X, filt=0.1, binarize=False):
    """Removes columbns with < filt non-zero entries"""
    if not binarize:
        return X.loc[:, X.sum(axis=0)/X.shape[0] >= filt]
    return X.loc[:, (X!=0).sum(axis=0)/X.shape[0] >= filt]

def count_filter(X, filt=1e-8):
    """ sets entries < filter to 0"""
    return (X >= filt) * X


def get_num_nodes(gfa_file):
    with open(gfa_file) as f:
        nodes = {r.nid for r in ps.parse_gfa(f) if r.type == ps.GFATypes.S}

    # make sure copangraph nodes are numbers 1:N
    for i in range(1, len(nodes)+1):
        if str(i) not in nodes:
            print(i)
    return len(nodes)


def calc_adjM(gfa):
    num_nodes = get_num_nodes(gfa)
    adjM = lil_matrix((num_nodes, num_nodes), dtype=bool)
    with open(gfa) as f:
        for e in ps.parse_gfa(f):
            if e.type != ps.GFATypes.L:
                continue
            adjM[int(e.l_nid)-1, int(e.r_nid)-1] = True
            adjM[int(e.r_nid)-1, int(e.l_nid)-1] = True
    return adjM
            

def compute_ct_mat(y, X):
    y = y.values[:, np.newaxis]
    X = X.values

    ct = np.array([
        [(y == 0) & (X == 0), (y == 0) & (X == 1)],
        [(y == 1) & (X == 0), (y == 1) & (X == 1)]
    ]).sum(axis=2)
    return ct

def compute_fisher_tests(y, X, test=stats.fisher_exact):
    results = dict()
    contingency_table_matrix = compute_ct_mat(y, X)
    print(contingency_table_matrix.shape)
    for i in tqdm.tqdm(range(contingency_table_matrix.shape[2])):
        ct = contingency_table_matrix[:,:,i]
        o,  p = test(ct)
        results[X.columns[i]] = {'statistic': o, 'pval':p}
    results = pd.DataFrame(results)
    return results.T

def compute_mwu_tests(y, X, test=stats.mannwhitneyu):

    # split the data into outcomeT (outcome==True) and outcomeF (outcome==False)
    outcomeT = X.loc[y.apply(lambda x: True if x == 1 else False), :]
    outcomeF = X.loc[y.apply(lambda x: False if x == 1 else True), :]
    mwu = test(outcomeT, outcomeF)
    res = pd.DataFrame([mwu.statistic, mwu.pvalue], index=['statistic', 'pval'], columns=X.columns).T
    return res


def fdr_correct(res, method='fdr_bh'):
    rejected, pval_corr, _, _ = multipletests(res.pval, method=method)
    res['rejected'] = rejected
    res['pval_corr'] = pval_corr
    return res
