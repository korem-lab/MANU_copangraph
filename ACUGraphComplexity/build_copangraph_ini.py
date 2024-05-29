import pandas as pd
COASM_TABLE = './data/ACUGraphComplexity/coassembly_table.csv'

INI_TMPLT = """
# enter all config in 'key = value' format
app_name = cpg
log_file_dir = logs
log_level = 0 
log_to = 0 
sample_list = /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/coasm_VARN_rep_VARREP.list
graph_name = coasm_VARN_rep_VARREP_ms_VARMS_dt_VARDT
out_dir = /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/
divergence_threshold = VARDT
num_threads = 64
max_separation = 150
window_size = 10
kmer_size = 15
min_homology_overlap = VARMS
min_contiguity_overlap = 60
max_jump = 200 
high_freq_kmer_filter = 1e-5 
fasta_file_ext = .fasta 
gfa_file_ext = .gfa
node_color_file_ext = .ncolor.csv
edge_color_file_ext = .ecolor.csv
extended_contigs = true
sensitive_mode = true
"""
if __name__ == '__main__':
    
    coassembly_table = pd.read_csv(COASM_TABLE)
    
    # construct ini and sample_list file, and coassembly pathnames
    groups = coassembly_table.groupby(by=['replicate', 'N'])
    for ((rep, n), df) in groups:
        rep = int(rep)
        n = int(n)
        print(rep, n)
        with open(f'./data/ACUGraphComplexity/coasm_{n}_rep_{rep}.list', 'w') as f:
            f.write('\n'.join(list(df.path)) + '\n')
        for ms in [100, 1000, 5000]:
            with open(f'./data/ACUGraphComplexity/coasm_{n}_rep_{rep}_ms_{ms}.ini', 'w') as f:
                out = INI_TMPLT.replace('VARN', str(n)).replace('VARMS', str(ms)).replace('VARREP',str(rep)).replace('VARDT', '0.02')
                f.write(out + '\n')
        for dt in ['0.005', '0.01', '0.02', '0.05', '0.1']:
            with open(f'./data/ACUGraphComplexity/coasm_{n}_rep_{rep}_ms_1000_dt_{dt}.ini', 'w') as f:
                out = INI_TMPLT.replace('VARN', str(n)).replace('VARMS', '1000').replace('VARREP',str(rep)).replace('VARDT', dt)
                f.write(out + '\n')


