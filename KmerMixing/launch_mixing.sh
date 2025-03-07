#!/bin/bash 
#SBATCH --job-name=kmer-mixing
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=m006
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/mixing-%x-%j.log
python run.py $KMER
