import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == '__main__':
    copan = pd.read_csv('./test_data/dummy_test_copan.csv')
    coasm = pd.read_csv('./test_data/dummy_test_coasm.csv')
    print(copan.shape)
    copan = copan.loc[ copan.sample_count > 1, :]
    print(copan.shape)
    print(coasm.shape)
    coasm= coasm.loc[coasm.sample_count > 1, :]
    print(coasm.shape)
    #sns.kdeplot(x=copan.sample_count, y=copan.node_count)
    cpg_coasm_prop = (copan.node_count < copan.sample_count).sum() / len(copan.node_count)
    coa_coasm_prop = (coasm.node_count < coasm.sample_count).sum() / len(coasm.node_count)
    print('copan mixing prop: ', cpg_coasm_prop)
    print('coasm mixing prop: ', coa_coasm_prop)
    plt.show()
    