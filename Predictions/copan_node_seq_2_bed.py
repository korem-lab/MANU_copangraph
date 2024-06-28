import sys
import os
import pandas as pd
import parse_seq
import glob
import re
DATA = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/mdro'

if __name__ == '__main__':
    
    gfas = glob.glob(os.path.join(DATA, '*.gfa'))  
    sns = pd.read_csv(os.path.join(DATA, 'mdro_pos.list'), header=None)
    for gfa in gfas:
        print(gfa, '...')
        u, v = None, None
        if '+' in gfa or '-' in gfa:
            u, v = re.findall('([0-9]+)[\\-\\+]_([0-9]+)[\\+\\-]\\.gfa', gfa)[0]
            print(u, v)
        else:
            u = re.findall('([0-9]+)\\.gfa', gfa)[0]
            print(u)
            
        pref = gfa.replace('.gfa', '')
        records = list()
        with open(gfa) as f:
            for s in parse_seq.parse_gfa(f):
                if s is None or s.type != parse_seq.GFATypes.S:
                    continue
                if s.nid == u or s.nid == v:
                    rest = s.rest[0].split(':')
                    cn, sn, lb, rb, ori = rest[2], rest[3], rest[5], rest[6], rest[7]
                    full_name = sns.loc[int(sn), 0]
                    records.append((f'{cn}_{sn}_{full_name}', lb, rb))

        with open(pref + '.bed', 'w') as f:
            for r in records:
                f.write('\t'.join(r) + '\n')

    