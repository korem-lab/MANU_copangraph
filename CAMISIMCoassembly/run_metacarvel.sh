#!/bin/bash 
#SBATCH --job-name=mcvl
#SBATCH --time=72:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --account pmg
#SBATCH --output=%x-%j.log
MCVL=/burg/pmg/users/ic2465/MetaCarvel/run.py

reads=$INPUT
DATA=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/CAMISIMCoassembly
bn=`basename ${reads%.txt}`
SCRATCH=/pmglocal/ic2465/${bn}_mcvl
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

cd $SCRATCH
pigz pooled_R1.fastq
pigz pooled_R2.fastq
bowtie2-build --threads 8 ./megahit/final.contigs.fa idx #build the index

bowtie2 -x idx -U pooled_R1.fastq.gz | samtools view -bS - | samtools sort - -o alignment_1.bam #align first reads
bowtie2 -x idx -U pooled_R2.fastq.gz | samtools view -bS - | samtools sort - -o alignment_2.bam #align second reads
samtools merge alignment_total.bam alignment_1.bam alignment_2.bam #merge the alignments 
samtools sort -n alignment_total.bam -o alignment.bam #sort by read names 
python $MCVL --assembly ./megahit/final.contigs.fa --mapping alignment.bam --dir mcvl

rm pooled_R1.fastq.gz
rm pooled_R2.fastq.gz
cd ..
mv $SCRATCH $DATA
