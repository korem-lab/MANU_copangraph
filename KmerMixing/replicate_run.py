import glob
import sys
import os
import mixing
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

GENOMEDIR='/manitou/pmg/projects/korem_lab/Projects/bacmeta_for_copan_v2/bacmeta_SIDNPOP25NBAC2'
MEGAHITDIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/replicates2/megahit'
COPANDIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/replicates2/copangraph'
OUTDIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/replicates2/'

if __name__ == '__main__':
    replicate = sys.argv[1]
    gen= sys.argv[2]

    # make kmer sets
    genome_list = glob.glob(GENOMEDIR + f'/{replicate}/fastas/Generation_{gen}' + f'/*.fasta')
    print(genome_list)
    kmer_sets = mixing.get_kmer_sets(genome_list)
    
    # make coasm
    seq_in_genome = list()
    for fg in glob.glob(MEGAHITDIR + f'/{replicate}_{gen}_pooled/intermediate_contigs/*.fastg'):
        print(fg)
        df_mix = mixing.seq_in_genome_mixing(
            fg, 
            genome_list,
            fg,
            kmer_sets, 
            candidate_its=25
        )
        seq_in_genome.append(df_mix)

    for fg in glob.glob(COPANDIR + f'/{replicate}_{gen}*.gfa'):
        print(fg)
        df_mix = mixing.seq_in_genome_mixing(
            fg, 
            genome_list,
            fg,
            kmer_sets, 
            candidate_its=25
        )
        seq_in_genome.append(df_mix)

    seq_in_genome = pd.concat(seq_in_genome)
    seq_in_genome.to_csv(f'{OUTDIR}/{replicate}_{gen}.csv')
