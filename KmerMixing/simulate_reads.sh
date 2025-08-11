#!/bin/bash 
#SBATCH --job-name=sim-reads
#SBATCH --time=0:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --exclude=m002,m006
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/KmerMixing/sim-reads-%x-%j.log

if [[ $FASTA == *"ZvJQaO"* ]]; then
	rep=ZvJQaO
elif [[ $FASTA ==  *"NKKBqt"* ]]; then
	rep=NKKBqt
elif [[ $FASTA == *"ECValN"* ]]; then
	rep=ECValN
else
	echo "NOT A REPLICATE"
fi

if [[ $FASTA == *"9000"* ]]; then
	gen=9000
elif [[ $FASTA == *"5000"* ]]; then
	gen=5000
elif [[ $FASTA == *"1000"* ]]; then
	gen=1000
else
	echo "NOT A GEN"
fi

OUTDIR=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/KmerMixing/bac_meta/replicates2/fastq/
base=`basename $FASTA`
art_illumina -ss HS25 -rs 42 -p -m 400 -s 10 -f 30 -i $FASTA -l 150 -o ${OUTDIR}/${rep}_${gen}_${base%.*}_

