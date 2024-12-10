import os
import sys
import glob
from subprocess import run
import parse_seq as ps
import pandas as pd
import numpy as np
import yaml
import shutil
from concurrent.futures import ProcessPoolExecutor
#moms-pi
#READ_DIR = '/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/HGF'
#OUT_DIR = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/order_invariance/moms-pi'
#CONTIG_DIR = '/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/EXT/megahit'
#GFA_FILE = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/order_invariance/moms-pi/moms-pi.gfa'
#GRAPH_SAMPLES = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/order_invariance/moms-pi/samples.csv'
#MAP_SAMPLES = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/order_invariance/moms-pi/samples.csv'
#READ_COUNTS = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/order_invariance/moms-pi/samples_read_counts.csv'

## acu
#READ_DIR = '/manitou/pmg/projects/korem_lab/Projects/ACU_PLT/2022_ACU_PLT/HGF2/'
#OUT_DIR = '/burg/pmg/users/ic2465/order_invariance/node_abund_scratch/'
#CONTIG_DIR = '/manitou/pmg/projects/korem_lab/Projects/ACU_PLT/2022_ACU_PLT/EXT/megahit/'
#GFA_FILE = '/burg/pmg/users/ic2465/order_invariance/persistence_sorted_qderived.gfa'
#GRAPH_SAMPLES = '/burg/pmg/users/ic2465/order_invariance/node_abund_scratch/graph_samples.csv'
#MAP_SAMPLES = '/burg/pmg/users/ic2465/order_invariance/node_abund_scratch/map_samples.csv'
#READ_COUNTS = '/manitou/pmg/projects/korem_lab/Projects/ACU_PLT/2022_ACU_PLT/HGF2/ReadCounts.csv'
#MIN_CONTIG_1000=1000
MIN_CONTIG_500=500
#MIN_CONTIG_250=250

def get_num_nodes(gfa_file):
    with open(gfa_file) as f:
        nodes = {r.nid for r in ps.parse_gfa(f) if r.type == ps.GFATypes.S}

    # make sure copangraph nodes are numbers 1:N
    for i in range(1, len(nodes)+1):
        if str(i) not in nodes:
            print(i)
    return len(nodes)


def process_sample(x):
        s, reads_dir, threads, out_dir = x
        s_1 = os.path.join(reads_dir, f'{s}_R1.fastq.gz')
        s_2 = os.path.join(reads_dir, f'{s}_R2.fastq.gz')
        assert os.path.isfile(s_1) and os.path.isfile(s_2)
        print('mapping', s_1, s_2)
        run(f'bowtie2 --threads {threads} -x {out_dir}/idx -1 {s_1} -2 {s_2} | samtools view -b -o {out_dir}/{s}_mapping.bam', shell=True)
        run(f'samtools sort -m 1G --threads {threads//2} -o {out_dir}/{s}_sorted_mapping.bam {out_dir}/{s}_mapping.bam', shell=True)
        run(f'samtools index --threads {threads} -o {out_dir}/{s}_sorted_mapping.bai {out_dir}/{s}_sorted_mapping.bam', shell=True)
        r = run(f'samtools idxstats {out_dir}/{s}_sorted_mapping.bam > {out_dir}/{s}_idxstats.csv', shell=True)
        return r.returncode

def map_sample(gfa_file, contigs_dir, reads_dir, graph_samples, map_samples, out_dir, num_reads_series, threads=1):
    """
     gfa_file: a copangraph made from samples in the graph_samples list
    """
    # Pool the contigs used to construct the graph into a single file
    print('Pooling contigs.')
    cntgs = [f'{contigs_dir}/{e}/final.contigs.fa' for e in graph_samples]
    assert all(os.path.isfile(e) for e in cntgs)
    cntg_out = f'{out_dir}/all.contigs.fa'
    fout = open(cntg_out, 'w')
    for s, cntg in zip(graph_samples, cntgs):
        with open(cntg) as f:
            for r in ps.parse(f, ps.Fasta):
                r.hdr = f'{s}:{r.hdr.split()[0]}'
                r.write(fout)
    fout.close()

    ## construct index
    print('Constructing index.')
    if (len(glob.glob(f'{out_dir}/idx*')) == 0):
        run(f'bowtie2-build --threads {threads*2} {cntg_out} {out_dir}/idx', shell=True)
    else:
        print("Index exists... will not recompute")

    # construct depth profiles for the test samples
    with ProcessPoolExecutor(max_workers=4) as ex:
        ex.map(process_sample, zip(map_samples, [reads_dir]*len(map_samples), [threads]*len(map_samples), [out_dir]*len(map_samples)))
        # map reads
    #for s in [e for e in map_samples if 'SRR6744867' == e]:
    #    s_1 = os.path.join(reads_dir, f'{s}_R1.fastq.gz')
    #    s_2 = os.path.join(reads_dir, f'{s}_R2.fastq.gz')
    #    print('mapping', s_1, s_2)
    #    run(f'bowtie2 --threads {threads} -x {out_dir}/idx -1 {s_1} -2 {s_2} | samtools view -b -o {out_dir}/{s}_mapping.bam', shell=True)
    #    run(f'samtools sort -m 1G --threads {threads//2} -o {out_dir}/{s}_sorted_mapping.bam {out_dir}/{s}_mapping.bam', shell=True)
    #    run(f'samtools index --threads {threads} -o {out_dir}/{s}_sorted_mapping.bai {out_dir}/{s}_sorted_mapping.bam', shell=True)
    #    run(f'samtools idxstats {out_dir}/{s}_sorted_mapping.bam > {out_dir}/{s}_idxstats.csv', shell=True)

    # make count matrix
    idx_dfs = [pd.read_csv(os.path.join(out_dir, f'{e}_idxstats.csv'), delimiter='\t', header=None, index_col=0).sort_index() for e in map_samples]
    for i in range(len(idx_dfs)-1):
        assert(all(idx_dfs[i].index == idx_dfs[i+1].index))
    cnts = pd.concat((e.loc[:, [2,3]].sum(axis=1) for e in idx_dfs), axis=1)
    cnts.columns = map_samples
    cnts['length'] = idx_dfs[0].loc[:, 1]  # get lengths
    idx_dfs = list() 
    print(cnts)

    # construct node features
    num_nodes = get_num_nodes(gfa_file)
    node_features = np.full((len(map_samples), num_nodes), 0.0)
    short_500= np.full((num_nodes,), False)
    with open(gfa_file) as f:
        for r in ps.parse_gfa(f):
            if r.type != ps.GFATypes.S:
                continue
            key = f'{graph_samples[r.sample_name]}:{r.contig_name}'
            clen = cnts.loc[key, 'length']
            if clen < MIN_CONTIG_500:
                short_500[int(r.nid)-1] = True 
            scale = abs(r.rb - r.lb) / clen
            node_features[:, int(r.nid)-1] += cnts.loc[key, map_samples].values * scale

    # library size normalization
    node_features = node_features / num_reads_series[map_samples].values[:, np.newaxis]
    node_features = pd.DataFrame(node_features, index=map_samples, columns=range(1, num_nodes+1))

    # remove nodes deriving from short contigs
    return node_features.T, node_features.loc[:, ~short_500].T,


if __name__ == '__main__':

    yml = sys.argv[1]
    assert(os.path.isfile(yml))
    with open(yml) as f:
        yml = yaml.safe_load(f)
    os.makedirs(yml['out_dir'], exist_ok=True)
    shutil.rmtree(yml['scratch_dir'], ignore_errors=True)
    os.makedirs(yml['scratch_dir'], exist_ok=True)
    graph_samples = pd.read_csv(yml['graph_samples'], header=None)[0]
    print('num_graph_sample: ',graph_samples.shape)
    map_samples = pd.read_csv(yml['map_samples'], header=None)[0]
    num_reads = pd.read_csv(yml['read_counts'], index_col=0, header=None)
    num_reads = num_reads.loc[map_samples, 1]
    print('map_samples: ', num_reads.shape)
    nf, nf_s500 = map_sample(yml['gfa_file'], yml['contig_dir'], yml['read_dir'], graph_samples, map_samples, yml['scratch_dir'], num_reads, threads=8)
    print('total node features', nf.shape)
    print('node_features_500', nf_s500.shape)
    nf.to_pickle(os.path.join(yml['out_dir'], f'X_full.T.pkl'))
    nf_s500.to_pickle(os.path.join(yml['out_dir'], f'X.T.pkl'))
    
