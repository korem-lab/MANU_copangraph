import os
import pandas as pd
import yaml

TASK='C'
EXT = '/manitou/pmg/projects/korem_lab/Projects/ACU_PLT/mmmbp2/tmp_path2/EXT/extended_contigs/'
SPLITS_FILE = f'/manitou/pmg/projects/korem_lab/Projects/MANU_copangraph/Predictions/mdro-graph-cv/task{TASK}/acu_123samp_longform_splits_nested.csv'
SEQ_DIV='0.02'
ASYMMETRIC='true'
FEATURE_DIR=f'/manitou/pmg/projects/korem_lab/Projects/MANU_copangraph//Predictions/mdro-graph-cv/ACU_TASK{TASK}_ASYM'
BASH_COPAN_TMPLT = f"""#!/bin/bash 
#SBATCH --job-name=cpn-XNM
#SBATCH --time=1:30:0
#SBATCH --mem=96G
#SBATCH --account=pmg
#SBATCH --cpus-per-task=4
#SBATCH --output={FEATURE_DIR}/%i-%j.log
/usr/bin/time -v /burg/pmg/users/ic2465/copangraph/bin/release/copangraph -s XSLIST -g graph -o XOUT -t 4 -d XSD -ms 250 -mj 750 -as XASYM
pigz XOUT/*.fasta XOUT/*.gfa
"""
READ_DIR ='/manitou/pmg/projects/korem_lab/Projects/ACU_PLT/mmmbp2/tmp_path2/UNT/'
CONTIG_DIR='/manitou/pmg/projects/korem_lab/Projects/ACU_PLT/mmmbp2/tmp_path2/EXT/megahit/'
READ_COUNTS='/manitou/pmg/projects/korem_lab/Projects/ACU_PLT/mmmbp2/df_path2/ReadCounts.csv'
TASK_SET=f'ACU_TASK{TASK}_ASYM'

BASH_MAP_TMPLT="""#!/bin/bash 
#SBATCH --job-name=abnd
#SBATCH --time=96:0:0
#SBATCH --mem=32G
#SBATCH --account=pmg
#SBATCH --cpus-per-task=4
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
        with open(os.path.join(key_dir, 'run_copan.sh'), 'w') as f:
            f.write(BASH_COPAN_TMPLT.replace('XNM', key[0]).replace('XOUT', key_dir).replace('XSLIST', sample_list).replace('XSD', SEQ_DIV).replace('XASYM', ASYMMETRIC))

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
            'read_counts': READ_COUNTS,
            'ext': '' 
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
            f.write(BASH_MAP_TMPLT.replace('XYML', yml_fl))
            



    
