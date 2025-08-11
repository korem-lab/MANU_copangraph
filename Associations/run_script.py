import sys
import utils
import pandas as pd

if __name__ == '__main__':

    abnd, outcome = sys.argv[1], sys.argv[2]
    #df = utils.run_adonis(abnd, outcome)
    abnd = pd.read_pickle(abnd).T
    abnd.index = abnd.index.astype(int)
    outcome = pd.read_csv(sys.argv[2], index_col=0)
    abnd = abnd.loc[outcome.index, :]
    print(outcome)
    print(abnd)
    print(abnd.shape)
    abnd = utils.count_filter(abnd, filt=10)
    abnd = utils.abundance_filter(abnd, filt=0.1, binarize=True)
    print(abnd.shape)
    if 'abx-load' in sys.argv[2]:
        res = utils.vectorized_spearman_with_p(outcome, abnd)
        res = utils.fdr_correct(res)
        print(res.loc[res.rejected & (res.pval_corr <= 0.01), :])
        res.to_csv('../data/Associations/abx-load/abx-load_02.spearman.csv')
    if 'ACU' in sys.argv[2]:
        outcome = outcome.outcome.apply(lambda x: 1 if x == 'Persistence' else 0)
        print(outcome)
        res = utils.compute_mwu_tests(outcome, abnd)
        res = utils.fdr_correct(res)
        print(res.loc[res.rejected, :])
        print(res.loc[res.pval < 0.001, :])
        res.to_csv(sys.argv[1].replace('pkl', 'mwu.csv'))

