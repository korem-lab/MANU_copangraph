import sys
import os

bs = """#!/bin/bash 
#SBATCH --job-name=bpeval
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --nodelist=m013
#SBATCH --account pmg
#SBATCH --output=/pmglocal/ic2465/%x-%j.log

python breakpoint.py X"""
file_name = os.path.basename(sys.argv[1])
file_name = os.path.splitext(file_name)[0]
with open(f'{file_name}.sh', 'w') as fo:
    fo.write(bs.replace('X', sys.argv[1]))
