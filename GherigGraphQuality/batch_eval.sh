#!/bin/bash 
#SBATCH --job-name=grg-co
#SBATCH --time=4:00:00
#SBATCH --mem=45G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --exclude=m002
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/GherigGraphQuality/log_%x-%j.log
/usr/bin/time -v python /burg/pmg/users/ic2465/Projects/MANU_copangraph/GherigGraphQuality/run.py $ANALYSIS

