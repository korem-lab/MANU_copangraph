#!/bin/bash 
#SBATCH --job-name=grig-co
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --exclude=m007
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/GherigGraphQuality/log_%x-%j.log
python /burg/pmg/users/ic2465/Projects/MANU_copangraph/GherigGraphQuality/run.py $ANALYSIS

