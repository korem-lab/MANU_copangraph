#!/bin/bash 
#SBATCH --job-name=ms
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/mdro/prokka-%x-%j.log
OUTDIR=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/mdro/prokka
DATA=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/Predictions/mdro/

base=`echo $FILE | rev | cut -b 7- | rev`
echo prokka --force --cpus 4 --outdir $OUTDIR --prefix $base $DATA/$FILE
prokka --force --cpus 4 --outdir $OUTDIR --prefix $base $DATA/$FILE