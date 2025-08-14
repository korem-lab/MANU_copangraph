import sys
import os
import glob
import pandas as pd
import Utils.parse_seq as ps

PAIRED_DATA='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/paired_data.csv'
MAG_DIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/drep_representatives_99'
SAMPLES_3 = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/3_sample_names.csv'
SAMPLES_6 = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/6_sample_names.csv'
SAMPLES_9 = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/9_sample_names.csv'
OUTDIR = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/subgraph_analysis' 


def construct_reference_coassembly(samples, pd_dict, coasm_sz):
    with open(samples) as f:
        data = [e.strip() for e in f]

    # get relevant mags
    print(data)
    mags = list()
    for sr in data:
        lr = pd_dict[sr]
        mags += glob.glob(os.path.join(MAG_DIR, f'{lr}*.fa'))
    print('num_mags', len(mags))

    with open(os.path.join(OUTDIR, f'{coasm_sz}_sample_reference.fasta'), 'w') as out:
        for mag in mags:
            with open(mag) as f:
                for fa in ps.parse(f, ps.Fasta):
                    fa.hdr = os.path.splitext(os.path.basename(mag))[0]
                    fa.write(out)

def construct_reference_sample(sample, pd_dict):
    lr = pd_dict[sample]
    mags = glob.glob(os.path.join(MAG_DIR, f'{lr}*.fa'))
    print(s, 'num_mags', len(mags))
    with open(os.path.join(OUTDIR, f'{lr}_reference.fasta'), 'w') as out:
        for mag in mags:
            with open(mag) as f:
                for fa in ps.parse(f, ps.Fasta):
                    fa.hdr = os.path.splitext(os.path.basename(mag))[0]
                    fa.write(out)
    

if __name__ == '__main__':

    paired_data = pd.read_csv(PAIRED_DATA)
    print(paired_data)
    pd_dict = dict()
    for i in paired_data.index:
        pd_dict[paired_data.loc[i,'short']] = paired_data.loc[i, 'long']
    
    # 3-sample
    #construct_reference_coassembly(SAMPLES_3, pd_dict, 3)
    #construct_reference_coassembly(SAMPLES_6, pd_dict, 6)
    #construct_reference_coassembly(SAMPLES_9, pd_dict, 9)

    with open(SAMPLES_9) as f:
        samples = [e.strip() for e in f]
        for s in samples:
            construct_reference_sample(s, pd_dict)
    
    
    