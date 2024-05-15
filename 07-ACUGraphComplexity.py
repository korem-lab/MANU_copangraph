import os 
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

DATA = './data/ACUGraphComplexity'
if __name__ == '__main__':
   
    df = pd.read_csv(os.path.join(DATA, 'complexity_table.csv'))
   
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
        nodes = df.loc[df.metric == m, :]
        sns.lineplot(x=nodes.n_samples, y=nodes.value, hue=nodes.assembler)
        plt.tight_layout()
        plt.savefig(os.path.join(DATA, f'{m}_complexity.pdf', dpi=1400, bbox_inches='tight')
   