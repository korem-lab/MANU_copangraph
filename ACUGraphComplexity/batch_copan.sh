#!/bin/bash 
#SBATCH --job-name=cpg
#SBATCH --time=2:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=64
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity/log_%x-%j.log
data=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/ACUGraphComplexity
write=/pmglocal/ic2465
COPANGRAPH=/burg/pmg/users/ic2465/copangraph/bin/release/copangraph

/usr/bin/time -v $COPANGRAPH $INI 
