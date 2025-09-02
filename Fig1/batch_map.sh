#!/bin/bash 
#SBATCH --job-name=bt2
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/PanMetagenomeViz/log_%x-%j.log
data=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/PanMetagenomeViz
write=/pmglocal/ic2465/pmgviz_${SAMPLE}

rm -rf ${write}
mkdir ${write}
cp ${data}/${SAMPLE}_R1.fq.gz ${data}/${SAMPLE}_R2.fq.gz ${write}
cp ${data}/${SAMPLE}.final.contigs.fa ${write}
bowtie2-build --threads 8 ${write}/${SAMPLE}.final.contigs.fa ${write}/idx
bowtie2 --threads 8 -x ${write}/idx -1 ${write}/${SAMPLE}_R1.fq.gz -2 ${write}/${SAMPLE}_R2.fq.gz | samtools view -bS -h - > ${write}/${SAMPLE}.unsorted.bam
samtools sort -@ 8 -o ${write}/${SAMPLE}.bam ${write}/${SAMPLE}.unsorted.bam
cp ${write}/${SAMPLE}.bam ${data}