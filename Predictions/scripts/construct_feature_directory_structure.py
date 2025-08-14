import os
import pandas as pd
import yaml

EXT = '/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/EXT/extended_contigs'
SPLITS_FILE = '/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/moms-pi-graph-cv/momspi_longform_splits_nested.csv'
SEQ_DIV='0.02'
FEATURE_DIR='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/moms-pi-graph-cv/FEATURES_NESTED'
COPAN_INI_TMPLT = """# enter all config in 'key = value' format
app_name = cpg
log_file_dir = logs # path relative to the executable
log_level = 0 # possible values are 1ebug=0(includes Info and Error), Info=1(includes Error), Error=2, None=3
log_to = 0 # possible values are Console=0, File=1, ConsoleAndFile=2
sample_list = XSAMPLELIST
graph_name = graph # name of graph prefix
out_dir = XKEYDIR
divergence_threshold = XSD
num_threads = 8
max_separation = 250
window_size = 10
kmer_size = 15
min_homology_overlap = 1000
min_contiguity_overlap = 60
max_jump = 750 # largest allowed gap within an alignment
high_freq_kmer_filter = 1e-5 # remove the top 1-x percentile kmers where x is input
fasta_file_ext = .fasta 
gfa_file_ext = .gfa
node_color_file_ext = .ncolor.csv
edge_color_file_ext = .ecolor.csv
extended_contigs = true
sensitive_mode = true
asymmetric = false"""
READ_DIR ='/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/HGF/'
CONTIG_DIR='/manitou/pmg/projects/korem_lab/Projects/MOMSPI_COPAN_PREDICTION/EXT/megahit'
READ_COUNTS='/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/order_invariance/moms-pi/samples_read_counts.csv'
TASK_SET='NUMOM_NESTED_CV'

BASH_TMPLT="""#!/bin/bash 
#SBATCH --job-name=abnd
#SBATCH --time=8:0:0
#SBATCH --mem=32G
#SBATCH --account=pmg
#SBATCH --cpus-per-task=64
#SBATCH --exclude=m013
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/Predictions/scripts/compute_abundance_features/%i-%j.log
python /burg/pmg/users/ic2465/Projects/MANU_copangraph/Predictions/scripts/compute_abundance_features/map_samples.py XYML
"""
if __name__ == '__main__':

    # load splits, organize by fold, and sort by outcome
    splts = pd.read_csv(SPLITS_FILE)
    splts =  splts.groupby(['key']).apply(lambda x: x.sort_values(['outcome', 'type', 'sample'], ascending=False)).reset_index(drop='level_0')
    for key, df in splts.groupby(['key']):
        print(key[0])
        # make directory
        key_dir = os.path.join(FEATURE_DIR, key[0])
        os.makedirs(key_dir, exist_ok=True)

        # write out sample list file and ini
        df_trn = df.loc[(df.type == 'train') | (df.type == 'outer_train')]
        sample_list = os.path.join(key_dir, 'sample.list')
        df_trn['sample'].apply(lambda x: os.path.join(EXT, f'{x}.sorted.pe_ext.fasta')).to_csv(sample_list, index=None, header=None)
        with open(os.path.join(key_dir, 'copan.ini'), 'w') as f:
            f.write(COPAN_INI_TMPLT.replace('XSD', SEQ_DIV).replace('XKEYDIR', key_dir).replace('XSAMPLELIST', sample_list))

        # write graph samples and map samples
        df_trn['sample'].to_csv(os.path.join(key_dir, 'graph_samples.csv'), index=None, header=None)
        df['sample'].to_csv(os.path.join(key_dir, 'map_samples.csv'), index=None, header=None)

        # write yaml
        yml = {
            'read_dir': READ_DIR, 
            'out_dir': key_dir,
            'scratch_dir': f'/pmglocal/ic2465/{TASK_SET}/{key[0]}',
            'contig_dir': CONTIG_DIR,
            'gfa_file': os.path.join(key_dir, 'graph.gfa'),
            'graph_samples': os.path.join(key_dir, 'graph_samples.csv'),
            'map_samples': os.path.join(key_dir, 'map_samples.csv'),
            'read_counts': READ_COUNTS
        }
        yml_fl = os.path.join(key_dir, 'feature_construction.yaml')
        with open(yml_fl, 'w') as f:
            yaml.dump(yml, f, default_flow_style=False)

        # write y full metadata. 
        df.to_csv(os.path.join(key_dir, 'md.csv'), index=None)
        
        # write y outcome and sample
        df.loc[:, ['sample', 'StudyID', 'outcome']].to_csv(os.path.join(key_dir, 'outcome.csv'), index=None)

        # write bash template to map samples
        with open(os.path.join(key_dir, 'map.sh'), 'w') as f:
            f.write(BASH_TMPLT.replace('XYML', yml_fl))
            



    
