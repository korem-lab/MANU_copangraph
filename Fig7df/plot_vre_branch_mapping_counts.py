import pandas as pd
import seaborn as sns
from scipy.stats import barnard_exact
from matplotlib import pyplot as plt
plt.figure(figsize=(3,3))
if __name__ == '__main__':
    stats = pd.DataFrame(columns=['Fig', 'statistic', 'pval'])
    records = [
            ['C-branch', 'Clearance', 2],
            ['C-branch', 'Persistance', 1],
            ['P-branch', 'Clearance',1],
            ['P-branch', 'Persistance',12]
        ]
    Fig7d_df = pd.DataFrame(
        records, columns=['Branch', 'MappedVRELabel', 'N_mapped']
    )
    sns.barplot(
        x=Fig7d_df.Branch, y=Fig7d_df.N_mapped, hue=Fig7d_df.MappedVRELabel
    )
    plt.tight_layout()
    plt.savefig('../data/Fig7df/Fig7d.png')
    res = barnard_exact(
        [ [Fig7d_df.loc[0,'N_mapped'], Fig7d_df.loc[2,'N_mapped']],
          [Fig7d_df.loc[1,'N_mapped'], Fig7d_df.loc[3,'N_mapped']]
        ],alternative='greater'
    )
    stats.loc[len(stats),:] = ['Fig7d', res.statistic, res.pvalue]

    records = [
                ['C-branch', 'Clearance', 4],
                ['C-branch', 'Persistance', 9],
                ['P-branch', 'Clearance', 0],
                ['P-branch', 'Persistance', 9],
        ]
    Fig7f_df = pd.DataFrame(
        records, columns=['Branch', 'MappedVRELabel', 'N_mapped']
    )
    sns.barplot(
        x=Fig7f_df.Branch, y=Fig7f_df.N_mapped, hue=Fig7f_df.MappedVRELabel
    )
    plt.tight_layout()
    plt.savefig('../data/Fig7df/Fig7f.png')
    res = barnard_exact(
        [ [Fig7f_df.loc[0,'N_mapped'], Fig7f_df.loc[2,'N_mapped']],
          [Fig7f_df.loc[1,'N_mapped'], Fig7f_df.loc[3,'N_mapped']]
        ],alternative='greater'
    )
    stats.loc[len(stats),:] = ['Fig7f', res.statistic, res.pvalue]
    stats.to_csv('../data/Fig7df/stats.csv')