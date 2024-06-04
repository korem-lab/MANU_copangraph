import os 
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

DATA = './data/ACUGraphComplexity'

LABELS = {
    '0.005': 'Copangraph (0.005)',
    '0.01': 'Copangraph (0.01)',
    '0.02': 'Copangraph (0.02)',
    '0.05': 'Copangraph (0.05)',
    '0.1': 'Copangraph (0.1)',
    '0.2': 'Megahit graph',
    }
if __name__ == '__main__':
   
    df = pd.read_csv(os.path.join(DATA, 'complexity_results.csv'), index_col=0)
   
    ## melt copangraph ms to single column
    #non_copan = df.loc[df.assembler != 'copangraph', :]
    #copan = df.loc[df.assembler == 'copangraph', :]
    #copan.columns = ['tmp'] + list(copan.columns[1:]) # rename assembler col
    #copan.loc[:, 'assembler'] = copan.index.to_series().apply(lambda x: copan.loc[x, 'tmp'] + '_' + str(copan.loc[x, 'ms']))
    #copan = copan.drop(columns=['tmp'])
    #df = pd.concat([copan, non_copan])
   
    ## plot nodes
    #for m in ['nodes', 'edges', 'N50', 'N90']:
    #    plt.clf()
    #    m_df = df.loc[df.metric == m, :]
    #    sns.lineplot(x=m_df.n_samples, y=m_df.value, hue=m_df.assembler)
    #    plt.tight_layout()
    #    plt.savefig(os.path.join(DATA, f'{m}_complexity.pdf'), dpi=1400, bbox_inches='tight')
    #    
    mh = df.loc[df.assembler == 'megahit', :]
    mh.assembler = mh.assembler.apply(lambda x: f'{x}_graph')
    df_dt = pd.read_csv(os.path.join(DATA, 'complexity_results_dt.csv'), index_col=0)
    df_dt = pd.concat([df_dt, mh], ignore_index=True)
    df_dt = df_dt.fillna(0.2)
    df_dt['labels'] = df_dt.index.to_series().apply(
        lambda i: df_dt.loc[i, 'assembler'] + ' ' + str(df_dt.loc[i, 'dt'])
    )
    
    print(df_dt)
    # plot
    for m in ['nodes', 'edges', 'N50', 'N90']:
        plt.clf()
        m_df = df_dt.loc[df_dt.metric == m, :]
        sns.lineplot(x=m_df.n_samples, y=m_df.value, hue=m_df.dt, palette='flare')
        plt.tight_layout()
        legend = plt.gca().get_legend()
        handles, labels = legend.legendHandles, legend.get_texts()
        plt.legend(handles=handles, labels=[LABELS[e.get_text()] for e in labels])
        plt.savefig(os.path.join(DATA, f'{m}_complexity_dt.pdf'), dpi=1400, bbox_inches='tight')
    
   