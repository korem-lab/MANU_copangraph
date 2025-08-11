#!/bin/bash 
#SBATCH --job-name=cami-gc
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMGraphQuality/camisim_reads/quality_%x-%j.log
/usr/bin/time -v python ./CAMISIMGraphQuality/run.py $ANALYSIS
