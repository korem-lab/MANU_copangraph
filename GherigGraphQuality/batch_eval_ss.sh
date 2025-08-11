#!/bin/bash 
#SBATCH --job-name=grg-ss
#SBATCH --time=1:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/GherigGraphQuality/grgss_%x-%j.log
/usr/bin/time -v python /burg/pmg/users/ic2465/Projects/MANU_copangraph/GherigGraphQuality/run.py $ANALYSIS

