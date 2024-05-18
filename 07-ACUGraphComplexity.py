import os 
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

DATA = './data/ACUGraphComplexity'
if __name__ == '__main__':
   
    df = pd.read_csv(os.path.join(DATA, 'complexity_results.csv'), index_col=0)
   
    # melt copangraph ms to single column
    non_copan = df.loc[df.assembler != 'copangraph', :]
    copan = df.loc[df.assembler == 'copangraph', :]
    copan.columns = ['tmp'] + list(copan.columns[1:]) # rename assembler col
    copan.loc[:, 'assembler'] = copan.index.to_series().apply(lambda x: copan.loc[x, 'tmp'] + '_' + str(copan.loc[x, 'ms']))
    copan = copan.drop(columns=['tmp'])
    df = pd.concat([copan, non_copan])
   
    # plot nodes
    for m in ['nodes', 'edges', 'N50', 'N90']:
        plt.clf()
        m_df = df.loc[df.metric == m, :]
        sns.lineplot(x=m_df.n_samples, y=m_df.value, hue=m_df.assembler)
        plt.tight_layout()
        plt.savefig(os.path.join(DATA, f'{m}_complexity.pdf'), dpi=1400, bbox_inches='tight')
   