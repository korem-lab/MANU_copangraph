#!/bin/bash 
#SBATCH --job-name=mcvl
#SBATCH --time=5-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --account pmg
#SBATCH --output=%x-%j.log
MCVL=/burg/pmg/users/ic2465/MetaCarvel/run.py
r1_reads=/manitou/pmg/users/ic2465/Projects/MANU_copangraph_2022/data/GherigGraphQuality/coassemblies/${INPUT}_pooled_R1.fastq.gz
r2_reads=/manitou/pmg/users/ic2465/Projects/MANU_copangraph_2022/data/GherigGraphQuality/coassemblies/${INPUT}_pooled_R2.fastq.gz

DATA=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality
SCRATCH=/pmglocal/ic2465/gherig_${INPUT}_mcvl
rm -rf $SCRATCH
mkdir -p $SCRATCH
cp $r1_reads $SCRATCH/pooled_R1.fastq.gz
cp $r2_reads $SCRATCH/pooled_R2.fastq.gz

megahit -t 8 -1 $SCRATCH/pooled_R1.fastq.gz -2 $SCRATCH/pooled_R2.fastq.gz -o $SCRATCH/megahit

cd $SCRATCH
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
