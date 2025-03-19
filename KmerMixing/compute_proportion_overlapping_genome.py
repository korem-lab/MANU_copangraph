import glob
import sys
import os
import mixing
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

GENOMEDIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/gherig/bins/drep_representatives_ANI/'
MEGAHITDIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/gherig/coassembly/megahit/intermediate_contigs/'
COPANDIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/gherig/copangraph/'
OUTDIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/gherig/'

if __name__ == '__main__':
    ani = sys.argv[1]
    threads = int(sys.argv[2])
    GENOMEDIR = GENOMEDIR.replace('ANI', ani)
    # make kmer sets
    genome_list = glob.glob(GENOMEDIR + f'*.fa')
    seq_in_genome = list()
    for fg in glob.glob(MEGAHITDIR + f'*.fastg'):
        print(fg)
        df_mix = mixing.covered_positions(
            fg, 
            genome_list,
            candidate_its=25,
            threads=threads
        )
        seq_in_genome.append(df_mix)
        print(df_mix)

    for fg in glob.glob(COPANDIR + f'/gherig_copangraph_sd_*.gfa'):
        print(fg)
        df_mix = mixing.covered_positions(
            fg, 
            genome_list,
            candidate_its=25,
            threads=threads
        )
        seq_in_genome.append(df_mix)
        print(df_mix)

    covered_pos= pd.concat(seq_in_genome)
    covered_pos.to_csv(f'{OUTDIR}/ani_{ani}_covered_positions.csv')
