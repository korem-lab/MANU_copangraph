import pandas as pd
import sys


TMPLT = """asm_artifact_gap: 150
max_within_clust_distance: 200
window_size: 175
run_description: '1_sample_coasm_cocpg_alllr'
reference: './data/GherigGraphQuality/1_sample_pooled_flye_asm.fasta'
out_dir: './data/GherigGraphQuality/coassemblies/'
key: '1_sample_coasm_cocpg_alllr'
dataset: '1_sample_coasm_cocpg_alllr'
ASMS:
  ASM_0:
    assembler: 'megahit_contigs'
    assembly_file: './data/GherigGraphQuality/coassemblies/1_sample_coasm/final.contigs.fa'
  ASM_1:
    assembler: 'megahit_graph'
    assembly_file: './data/GherigGraphQuality/coassemblies/1_sample_coasm/k141.fastg'
  ASM_2:
    assembler: 'copangraph'
    assembly_file: './data/GherigGraphQuality/coassemblies/1_sample_cpg.gfa'
"""

TMPLT2 = """asm_artifact_gap: 150
max_within_clust_distance: 200
window_size: 175
run_description: '1_sample_ssasm_sscpg_sslr_SRXXX'
reference: './data/GherigGraphQuality/ASMXXX_flye_asm.fasta'
out_dir: './data/GherigGraphQuality/coassemblies/'
key: '1_sample_ssasm_sscpg_sslr_SRXXX'
dataset: '1_sample_SRXXX'
long_reads: 'LRXXX'
ASMS:
  ASM_0:
    assembler: 'megahit_contigs'
    assembly_file: './data/GherigGraphQuality/megahit/GHERIG_SRXXX/final.contigs.fa'
  ASM_1:
    assembler: 'megahit_graph'
    assembly_file: './data/GherigGraphQuality/megahit/GHERIG_SRXXX/k141.fastg'
  ASM_2:
    assembler: 'copangraph'
    assembly_file: './data/GherigGraphQuality/extractions/1_sample_SRXXX_cpg.gfa'
"""

def map_short_to_long(table, short):
    sample = table.loc[table.run == short, 'sample'].iloc[0]  # get sample name
    return table.loc[((table == sample)|(table == 'PACBIO_SMRT')).sum(axis=1) > 1, 'run'].iloc[0]

if __name__ == '__main__':
    
    if len(sys.argv) != 5:
        print('Usage: <exe> <long-short pair map> sample_list <num_samples_in_coasm> <mode>')
        sys.exit()
    mode = sys.argv[4]
    assert mode in [
        'ssasm_sscpg_sslr', 
        'coasm_cocpg_alllr',
    ]
    ls_map = sys.argv[1] 
    ls_map = pd.read_csv(ls_map, index_col=0)
    samples = sys.argv[2]
    samples = pd.read_csv(samples, header=None).loc[:,0]
    coasm_sz = sys.argv[3] 
    
    if mode == 'coasm_cocpg_alllr':
        with open(f'{coasm_sz}_coasm_cocpg_alllr.yaml', 'w') as f:
             out = TMPLT.replace('1_sample', f'{coasm_sz}_sample')
             print(out)
             f.write(out + '\n')
    elif mode == 'ssasm_sscpg_sslr':
        print(ls_map)
        for s in samples:
            s = s.replace('GHERIG_', '')
            l = map_short_to_long(ls_map, s)
            print(s, l)
            with open(f'{coasm_sz}_sample_ssasm_sscpg_sslr_{l}_{s}.yaml', 'w') as f:
                l_asm_pref = 'GHERIG_' + l
                out = TMPLT2.replace('ASMXXX', l_asm_pref).replace('SRXXX', s)\
                    .replace('LRXXX', l).replace('1_sample', f'{coasm_sz}_sample')
                print(out)
                f.write(out + '\n')
        
