import glob
import sys
import os
import mixing
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
from collections import defaultdict

if __name__ == '__main__':
    kmer = sys.argv[1]
    paths = defaultdict(defaultdict)
    for ANI, GEN in [('ANI90', '19000'), ('ANI95', '9000'), ('ANI97', '5000'), ('ANI99', '1000')]:
        paths[ANI]['FASTA'] = f'/manitou/pmg/projects/korem_lab/Projects/bacmeta_for_copan/bacmeta_SIDNPOP25NBAC2/CYbvDD/fastas/Generation_{GEN}'
        paths[ANI]['FASTQ'] =f'/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/{ANI}/fastq'
        paths[ANI]['MEGAHIT'] = f'/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/{ANI}/megahit'
        paths[ANI]['COPAN'] = f'/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/{ANI}/copan'
    
    
    # Run 'seq in genome' mixing
    seq_in_genome = list()
    kmer_sets = mixing.get_kmer_sets(glob.glob(paths['ANI95']['FASTA'] + '/*.fasta'))
    for fg in glob.glob(f'{paths["ANI95"]["MEGAHIT"]}/coasm/intermediate_contigs/*.fastg'):
        if kmer not in fg:
            continue
        print(kmer, fg)
        df_mix = mixing.seq_in_genome_mixing(
            fg, 
            glob.glob(paths['ANI95']['FASTA'] + '/*.fasta'), 
            fg,
            kmer_sets, 
            candidate_its=25
        )
        seq_in_genome.append(df_mix)
    for fg in glob.glob(f'{paths["ANI95"]["COPAN"]}/*.gfa'):
        if kmer not in fg:
            continue
        print(fg)
        df_mix = mixing.seq_in_genome_mixing(
            fg, 
            glob.glob(paths['ANI95']['FASTA'] + '/*.fasta'), 
            fg,
            kmer_sets, 
            candidate_its=25
        )
        seq_in_genome.append(df_mix)
    seq_in_genome = pd.concat(seq_in_genome)
    seq_in_genome.to_csv(f'seq_in_genome_{kmer}.csv')
