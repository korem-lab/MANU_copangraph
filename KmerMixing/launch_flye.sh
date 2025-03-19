#!/bin/bash 
#SBATCH --job-name=checkM
#SBATCH --time=4:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --exclude=m005
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/flye-%x-%j.log

flye --meta --threads 16 --pacbio-hifi ${SAMPLE} -o ${OUTDIR}

