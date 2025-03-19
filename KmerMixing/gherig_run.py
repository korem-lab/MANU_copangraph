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
    GENOMEDIR = GENOMEDIR.replace('ANI', ani)
    # make kmer sets
    genome_list = glob.glob(GENOMEDIR + f'*.fa')
    print(genome_list)
    kmer_sets = mixing.get_kmer_sets(genome_list)
    
    # make coasm
    seq_in_genome = list()
    for fg in glob.glob(MEGAHITDIR + f'*.fastg'):
        print(fg)
        df_mix = mixing.seq_in_genome_mixing(
            fg, 
            genome_list,
            fg,
            kmer_sets, 
            candidate_its=10
        )
        seq_in_genome.append(df_mix)

    for fg in glob.glob(COPANDIR + f'/gherig_copangraph_sd_*.gfa'):
        print(fg)
        df_mix = mixing.seq_in_genome_mixing(
            fg, 
            genome_list,
            fg,
            kmer_sets, 
            candidate_its=10
        )
        seq_in_genome.append(df_mix)

    seq_in_genome = pd.concat(seq_in_genome)
    seq_in_genome.to_csv(f'{OUTDIR}/ani_{ani}_mixing_results.csv')
