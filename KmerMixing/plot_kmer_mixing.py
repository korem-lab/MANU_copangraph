import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
DATA = './test_data/coasm_5_rep_0_'
if __name__ == '__main__':
    copan_1000_02 = pd.read_csv(DATA + 'copan_ms_1000_dt_02_kmix.csv')
    copan_100_02 = pd.read_csv(DATA + 'copan_ms_100_dt_02_kmix.csv')
    copan_1000_05 = pd.read_csv(DATA + 'copan_ms_1000_dt_05_kmix.csv')
    coasm_k141 = pd.read_csv(DATA + 'k141_kmix.csv')
    coasm_k59 = pd.read_csv(DATA + 'k59_kmix.csv')
    
    copan_1000_02 = copan_1000_02.loc[copan_1000_02.sample_count > 1, :]
    copan_100_02 = copan_100_02.loc[copan_100_02.sample_count > 1, :]
    copan_1000_05 = copan_1000_05.loc[copan_1000_05.sample_count > 1, :]
    coasm_k141 = coasm_k141.loc[coasm_k141.sample_count > 1, :]
    coasm_k59 = coasm_k59.loc[coasm_k59.sample_count > 1, :]
    
    
    #sns.kdeplot(x=copan.sample_count, y=copan.node_count)
    copan_1000_02_prop = (copan_1000_02.node_count < copan_1000_02.sample_count).sum() / len(copan_1000_02.node_count)
    copan_100_02_prop = (copan_100_02.node_count < copan_100_02.sample_count).sum() / len(copan_100_02.node_count)
    copan_1000_05_prop = (copan_1000_05.node_count < copan_1000_05.sample_count).sum() / len(copan_1000_05.node_count)
    coasm_k141_prop = (coasm_k141.node_count < coasm_k141.sample_count).sum() / len(coasm_k141.node_count)
    coasm_k59_prop = (coasm_k59.node_count < coasm_k59.sample_count).sum() / len(coasm_k59.node_count)
    
    print('copan_1000_02: ', copan_1000_02_prop)
    print('copan_100_02: ', copan_100_02_prop)
    print('copan_1000_05: ', copan_1000_05_prop)
    print('coasm_k141: ', coasm_k141_prop)
    print('coasm_k59: ', coasm_k59_prop)
    