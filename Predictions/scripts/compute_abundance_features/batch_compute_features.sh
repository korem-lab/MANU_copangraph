#!/bin/bash 
#SBATCH --job-name=momspi
#SBATCH --time=12:0:0
#SBATCH --mem=32G
#SBATCH --account=pmg
#SBATCH --cpus-per-task=4
#SBATCH --exclude=m013
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/Predictions/scripts/compute_abundance_features/log.txt
python map_samples.py
