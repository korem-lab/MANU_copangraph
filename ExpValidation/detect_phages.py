import sys
import os
import glob
import subprocess
import pandas as pd
from BashWriter import *

def detect_phages(indir,out_dir,batch_size,exclude_nodes):
    subprocess.run(f'mkdir {out_dir}/vibrant/',shell=True)
    subprocess.run(f'mkdir {out_dir}/genomad/',shell=True)
    bw = BashWriter(mem=32, hours=8, cpus=8,exclude=exclude_nodes)
    for file in glob.glob(os.path.join(indir, '*.fasta')):
        sample = os.path.basename(file).split('.')[0]
        cmd = f"""
        echo "-------------------------------------------------RUNNING VIBRANT-------------------------------------------------"
        source /burg/pmg/users/korem_lab/miniforge3/bin/activate VIBRANT
        cd $PMGLOCAL$
        cp {file} $PMGLOCAL$
        python /burg/pmg/users/korem_lab/Databases/vibrant/VIBRANT_run.py -i {sample}.fasta -no_plot -t 8
        mv VIBRANT* {out_dir}/vibrant/
        echo "-------------------------------------------------RUNNING GENOMAD-------------------------------------------------"
        source /burg/pmg/users/korem_lab/miniforge3/bin/activate genomad
        genomad end-to-end --cleanup --splits 8 {sample}.fasta {sample}_genomad /burg/pmg/users/korem_lab/Databases/genomad_db/
        mv {sample}_genomad/ {out_dir}/genomad/
        """
        bw.add(cmd)
    bw.run_batch(cmds_per_job=batch_size, name='subgraphs')


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python detect_phages.py <assembly_directory> <out_directory> <batch_size> <exclude_nodes>")
        sys.exit(1)
    
    ass_dir = sys.argv[1]
    out_dir = sys.argv[2]
    batch_size = int(sys.argv[3])
    exclude = sys.argv[4]
    detect_phages(ass_dir,out_dir,batch_size,exclude)
