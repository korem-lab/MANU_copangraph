import Utils.parse_seq as ps
import sys
import os
import glob
import re
import pandas as pd

def map_short_to_long(table, short):
    sample = table.loc[table.run == short, 'sample'].iloc[0]  # get sample name
    return table.loc[((table == sample)|(table == 'PACBIO_SMRT')).sum(axis=1) > 1, 'run'].iloc[0]

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: <exe> <long_short_read_map> <samples>', )
        sys.exit()
        
    ls_map = sys.argv[1] 
    ls_map = pd.read_csv(ls_map, index_col=0)
    samples = sys.argv[2]
    samples = pd.read_csv(samples, header=None).loc[:,0]
    coasm_sz = len(samples)
    asms = list()
    for s in samples:
        s = s.replace('GHERIG_', '')
        l = map_short_to_long(ls_map, s)
        asms.append(f'./data/GherigGraphQuality/GHERIG_{l}_flye_asm.fasta')
    with open(f'./data/GherigGraphQuality/{coasm_sz}_sample_pooled_flye_asm.fasta', 'w') as fout:
        for f in asms:
            srr = re.findall('GHERIG_(SRR[0-9]+)_', f)[0]
            with open(f) as fin:
                for fa in ps.parse(fin, ps.Fasta):
                    fa.hdr += '_' + srr
                    fa.write(fout)
    