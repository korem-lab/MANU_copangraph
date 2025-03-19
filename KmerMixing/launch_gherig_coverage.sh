#!/bin/bash 
#SBATCH --job-name=ovlp
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --exclude=m005
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/gherig-ovlp-%x-%j.log
python /burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/compute_proportion_overlapping_genome.py 99 32
