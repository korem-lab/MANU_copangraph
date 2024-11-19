#!/bin/bash 
#SBATCH --job-name=megahit
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --account pmg
#SBATCH --output=%x-%j.log


reads=$INPUT
DATA=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMCoassembly
bn=`basename ${reads%.txt}`
SCRATCH=/pmglocal/ic2465/$bn
rm -rf $SCRATCH
mkdir -p $SCRATCH
cat $reads | while read f; do
	if [[ "$f" == *"R1"* ]]; then
		zcat $f >> $SCRATCH/pooled_R1.fastq
	else
		zcat $f >> $SCRATCH/pooled_R2.fastq
	fi
done
megahit -t 8 -1 $SCRATCH/pooled_R1.fastq -2 $SCRATCH/pooled_R2.fastq -o $SCRATCH/megahit
rm $SCRATCH/pooled_R1.fastq
rm $SCRATCH/pooled_R2.fastq
mv $SCRATCH $DATA
