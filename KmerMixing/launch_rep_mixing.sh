#!/bin/bash 
#SBATCH --job-name=kmer-mixing
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --exclude=m005
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/repmix-%x-%j.log
python replicate_run.py ${r} ${g}
