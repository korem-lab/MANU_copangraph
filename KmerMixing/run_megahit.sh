#!/bin/bash 
#SBATCH --job-name=megahit
#SBATCH --time=2:0:0
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --exclude=m005
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/mh-%x-%j.log
FASTQ=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/replicates2/fastq
OUTDIR=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/replicates2/megahit/
megahit -t 8 -1 ${FASTQ}/${SAMPLE}_1.fastq -2 ${FASTQ}/${SAMPLE}_2.fastq -o ${OUTDIR}/${SAMPLE}

