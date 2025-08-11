from scipy import stats
import os
import numpy as np
import seaborn as sns
import pandas as pd
import parse_seq as ps
import tqdm
from subprocess import run
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

def compute_spearman_test(y, X, test=stats.spearmanr):
    print(X.values[:,0].shape)
    print(y.values.shape)
    stats, pvalues = list(), list()
    for i in tqdm.tqdm(range(X.shape[1])):
        s, p = test(X.values[:, i], y.values)
        stats.append(s)
        pvalues.append(p)
    res = pd.DataFrame([stats, pvalues], index=['statistic', 'pval'], columns=X.columns).T
    return res


def compute_mwu_tests(y, X, test=stats.mannwhitneyu, alternative='two-sided'):

    # split the data into outcomeT (outcome==True) and outcomeF (outcome==False)
    outcomeT = X.loc[y.apply(lambda x: True if x == 1 else False), :]
    outcomeF = X.loc[y.apply(lambda x: False if x == 1 else True), :]
    mwu = test(outcomeT, outcomeF, alternative=alternative)
    res = pd.DataFrame([mwu.statistic, mwu.pvalue], index=['statistic', 'pval'], columns=X.columns).T
    return res


def fdr_correct(res, method='fdr_bh'):
    rejected, pval_corr, _, _ = multipletests(res.pval, method=method)
    res['rejected'] = rejected
    res['pval_corr'] = pval_corr
    return res

def run_adonis(X_fl, outcome_fl):
    X_fl = X_fl.replace('.pkl', '.csv')
    if not os.path.exists(X_fl):
        print('making csv')
        x = pd.read_pickle(X_fl.replace('csv', 'pkl'))
        x.to_csv(X_fl)
    else:
        print('csv exists')
    assert(os.path.exists(X_fl))
    cmd = f'Rscript adonis.R {X_fl} {outcome_fl}'
    run(cmd, shell=True)
    adonis_fl = X_fl.replace('.csv', '_adonis.csv')
    adonis = pd.read_csv(adonis_fl)
    return adonis 

def run_lasso(X, y, split=False, penalty='l1', C=0.1):
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)
    if split:
        X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.3, random_state=101)
    else:
        X_train, X_test, y_train, y_test = X_scaled, X_scaled, y, y
    if penalty == 'elasticnet':
        solver = 'saga'
        l1_ratio = 0.5
    else:
        solver = 'liblinear'
        l1_ratio = None
    lasso = LogisticRegression(penalty=penalty, solver=solver, C=C, random_state=101, l1_ratio=l1_ratio)
    lasso.fit(X_train, y_train)
    y_pred_prob = lasso.predict_proba(X_test)[:, 1]
    auroc = roc_auc_score(y_test, y_pred_prob)
    print(f'AUROC: {auroc}')
    RocCurveDisplay.from_estimator(lasso, X_test, y_test)


def vectorized_spearman_with_p(y, X):
    # Rank the data
    X_ranked = np.apply_along_axis(stats.rankdata, 0, X)  # Rank columns of X
    y_ranked = stats.rankdata(y)  # Rank y

    # Center the ranks
    X_centered = X_ranked - np.mean(X_ranked, axis=0)
    y_centered = y_ranked - np.mean(y_ranked)

    # Compute Spearman correlations
    numerator = np.dot(X_centered.T, y_centered)
    denominator = np.sqrt(np.sum(X_centered**2, axis=0) * np.sum(y_centered**2))
    spearman_corrs = numerator / denominator

    # Compute p-values
    n = X.shape[0]
    t_stat = spearman_corrs * np.sqrt((n - 2) / (1 - spearman_corrs**2))
    p_values = 2 * stats.t.sf(np.abs(t_stat), df=n - 2)  # Two-tailed test
    res = pd.DataFrame([spearman_corrs, p_values], index=['statistic', 'pval'], columns=X.columns).T
    return res
