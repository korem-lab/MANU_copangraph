""" To be run on each split directory"""
import os
import sys
import pandas as pd
from glob import glob

DIR= '/burg/pmg/users/ic2465/Projects...'
FILE = 'all_dom_abund_pruned.tsv'

def construct_feature_table(files):
    features = list()
    for f in files:
        print(f)
        g=gen()
        data = pd.read_csv(f, delimiter='\t')
        data['name'] = f.split('/')[-2][:15]
        data['name'] = data['name'].apply(lambda x: f'{x}_{next(g)}')
        features.append(data)
    ret = pd.concat(features, ignore_index=True)
    return ret.T


if __name__ == '__main__':
    
    files = glob(os.path.join(DIR, 'GC*abund', FILE))
    for file in files:
        genome_name = f.split('/')[-2][:15]
        pd.read_csv(f, 
    features = construct_feature_table(files)
    features.to_pickle(f'momspi_spacegraphcats_features.pkl')
#        features.to_csv(f'{s}_spacegraphcats_features.csv', index=None)
        

