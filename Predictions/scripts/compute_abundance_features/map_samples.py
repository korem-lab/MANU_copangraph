import os
import sys
from subprocess import run
import parse_seq as ps
import pandas as pd
import numpy as np
TST_READS_DIR = '/home/izaak/Dropbox/Documents/Code/map_samples/mapped'
OUT_DIR = '/home/izaak/Dropbox/Documents/Code/map_samples/out_dir'
CONTIG_DIR = '/home/izaak/Dropbox/Documents/Code/map_samples/megahit'
GFA_FILE = '/home/izaak/Dropbox/Documents/Code/map_samples/G.gfa'
GRAPH_SAMPLES = ['30110', '30113', '30118']
MAP_SAMPLES = ['30119', '30192', '30193']

def get_num_nodes(gfa_file):
    with open(gfa_file) as f:
        nodes = {r.nid for r in ps.parse_gfa(f) if r.type == ps.GFATypes.S}

    # make sure copangraph nodes are numbers 1:N
    for i in range(1, len(nodes)+1):
        if str(i) not in nodes:
            print(i)
    return len(nodes)
        

def map_sample(gfa_file, contigs_dir, tst_reads_dir, s_trn, s_tst, out_dir, num_reads_series, threads=1):
    """
     gfa_file: a copangraph made from samples in the s_trn list
    """
    
    # Pool the contigs used to construct the graph into a single file
    print('Pooling contigs.')
    cntgs = [f'{contigs_dir}/{e}/final.contigs.fa' for e in s_trn]
    cntg_out = f'{out_dir}/all.contigs.fa'
    fout = open(cntg_out, 'w')
    for s, cntg in zip(s_trn, cntgs):
        with open(cntg) as f:
            for r in ps.parse(f, ps.Fasta):
                r.hdr = f'{s}:{r.hdr.split()[0]}'
                r.write(fout)
    fout.close()

    # construct index
    print('Constructing index.')
    run(f'bowtie2-build --threads {threads} {cntg_out} {out_dir}/idx', shell=True)

    # construct depth profiles for the test samples
    for s in s_tst:
        # map reads
        s_1 = os.path.join(tst_reads_dir, f'{s}_1.fastq.gz')
        s_2 = os.path.join(tst_reads_dir, f'{s}_2.fastq.gz')
        print('mapping', s_1, s_2)
        run(f'bowtie2 --threads {threads} -x {out_dir}/idx -1 {s_1} -2 {s_2} | samtools view -b -o {out_dir}/{s}_mapping.bam', shell=True)
        run(f'samtools sort --threads {threads} -o {out_dir}/{s}_sorted_mapping.bam {out_dir}/{s}_mapping.bam', shell=True)
        run(f'samtools index --threads {threads} -o {out_dir}/{s}_sorted_mapping.bai {out_dir}/{s}_sorted_mapping.bam', shell=True)
        run(f'samtools idxstats {out_dir}/{s}_sorted_mapping.bam > {out_dir}/{s}_idxstats.csv', shell=True)

    # make count matrix
    idx_dfs = [pd.read_csv(os.path.join(out_dir, f'{e}_idxstats.csv'), delimiter='\t', header=None, index_col=0).sort_index() for e in s_tst]
    for i in range(len(idx_dfs)-1):
        assert(all(idx_dfs[i].index == idx_dfs[i+1].index))
    cnts = pd.concat((e.loc[:, [2,3]].sum(axis=1) for e in idx_dfs), axis=1)
    cnts.columns = s_tst
    cnts['length'] = idx_dfs[0].loc[:, 1]  # get lengths
    idx_dfs = list() 
    print(cnts)

    # construct node features
    num_nodes = get_num_nodes(gfa_file)
    node_features = np.full((len(s_tst), num_nodes), 0.0)
    with open(gfa_file) as f:
        for r in ps.parse_gfa(f):
            if r.type != ps.GFATypes.S:
                continue
            key = f'{s_trn[r.sample_name]}:{r.contig_name}'
            clen = cnts.loc[key, 'length']
            scale = abs(r.rb - r.lb) / clen
            node_features[:, int(r.nid)-1] += cnts.loc[key, s_tst].values * scale

    # library size normalization
    node_features = node_features / cnts[s_tst].sum(axis=0).values[:, np.newaxis]
    return node_features


if __name__ == '__main__':

    num_reads = pd.Series([7379250, 9726297, 11807259], index=S_TST)
    print(map_sample(GFA_FILE, CONTIG_DIR, TST_READS_DIR, S_TRN, S_TST, OUT_DIR, num_reads, threads=8))
    
