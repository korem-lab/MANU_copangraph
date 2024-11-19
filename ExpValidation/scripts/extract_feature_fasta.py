import sys
import os
import re
import pandas as pd
import parse_seq as ps

if __name__ == '__main__':

    _, feature, fa, ncol, ecol, samples, out = sys.argv
    print(feature)
    # match edge
    match = re.findall('([0-9]+)([-\+])_([0-9]+)([-\+])', feature)
    print(match)
    if match:
        u, u_ori, v, v_ori, is_node = *match[0], False
    # match node
    else:
        u, v, is_node = re.findall('([0-9]+)', feature)[0], None, True

    print(feature, is_node)
    # extract fasta
    with open(fa) as fin, open(os.path.join(out, f'{feature}.fasta'), 'w') as fout:
        feature_fa = [e for e in ps.parse(fin, ps.Fasta) if (e.hdr.startswith(f'{u}:') or e.hdr.startswith(f'{v}:'))]
        print(len(feature_fa))
        for e in feature_fa:
            e.write(fout)
    # extract feature
    if is_node:
        df = pd.read_csv(ncol, index_col=0)
        feature_occ = df.loc[int(u), :]
    else:
        df = pd.read_csv(ecol, index_col=0)
        feature_occ = df.loc[f'{u}({u_ori}) -> {v}({v_ori})', :]
    samples = pd.read_csv(samples, header=None)
    feature_occ.index = samples.loc[:, 0]
    feature_occ.index.name = 'samples'
    
    # extract seq lengths
    lens = pd.Series([len(e.seq) for e in feature_fa], index=[e.hdr for e in feature_fa])

    # write
    feature_occ.to_csv(os.path.join(out, f'{feature}.occ.csv'), header=None)
    lens.to_csv(os.path.join(out, f'{feature}.seqlens.csv'), header=None)
