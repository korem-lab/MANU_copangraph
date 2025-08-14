"""
Copangraph should transform all input contigs into graph. No sequences should
be lost. This script quantifies lost base pairs. 
"""

import os 
import sys
import parse_seq as ps
import glob
from collections import defaultdict
import numpy as np
import tqdm

#GHERIG_SRR15489027
#GHERIG_SRR15489028
#GHERIG_SRR15489030
#GHERIG_SRR15489026
#GHERIG_SRR15489024
#GHERIG_SRR15489033
#GHERIG_SRR15489029
#GHERIG_SRR15489022
#GHERIG_SRR15489031

if __name__ == '__main__':
    contig_dir = '/manitou/pmg/projects/korem_lab/Projects/MANU_copangraph/data/GherigGraphQuality/old_EXT/megahit/GHERIG_SRR15489022/'
    gfa = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/oldsubgraph/9_sample_ssasm_sscpg_sslr_SRR15489022.gfa'
    fastafiles = glob.glob(os.path.join(contig_dir, 'final.contigs.fa'))
    print('Number of fasta files: ', len(fastafiles))


    # for each contig, build an zero array
    covered_positions = dict()
    for fastafile in fastafiles:
        print(fastafile)
        sample = fastafile.split('/')[-2]
        with open(fastafile) as f:
            for fa in ps.parse(f, ps.MegaHITFasta):
                tox = f'{sample}_{fa.cn}'
                covered_positions[tox] = np.zeros(len(fa.seq), dtype=np.uint8)

    # coverage vector calc
    with open(gfa) as f:
        for e in ps.parse_gfa(f):
            if e.type != ps.GFATypes.S:
                continue
            
            rest = e.rest.split(':')
            #tox = f'{rest[3]}_{rest[2]}'
            tox = f'GHERIG_SRR15489022_{rest[2]}'
            left, right = int(rest[5]), int(rest[6])
            covered_positions[tox][left: right] = 1

    # quantify coverage
    total_bp = 0
    total_covered_bp = 0
    total_contigs_with_missing_data = 0
    average_missing_bp = 0
    average_prop_missing = 0
    for k, v in tqdm.tqdm(covered_positions.items()):
        total_bp += len(v)
        total_covered_bp += v.sum()
        total_contigs_with_missing_data += 1 if (v == 0).any() else 0
        average_missing_bp += (v==0).sum()
        average_prop_missing += (v==0).sum() / float(len(v))
    
    average_missing_bp /= float(len(covered_positions))
    average_prop_missing /= float(len(covered_positions))
    print('Total bp: ', total_bp)
    print('Total covered bp: ', total_covered_bp)
    print('Proportion of total covered: ', float(total_covered_bp) / total_bp)
    print('=======')
    print('Total contigs: ', len(covered_positions))
    print('Total contigs missing data: ', total_contigs_with_missing_data)
    print('Proportion of contigs with missing data: ', total_contigs_with_missing_data / float(len(covered_positions)))
    print('=======')
    print('average_missing_bp: ', average_missing_bp)
    print('average_prop_missing: ', average_prop_missing)
