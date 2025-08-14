#!/bin/bash 
#SBATCH --job-name=mcvl
#SBATCH --time=6:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/ss-mcvl_%x-%j.log
OUTDIR=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/SSMCVL

TMPDIR=/pmglocal/tk2829/${SAMPLE}
rm -rf $TMPDIR
mkdir $TMPDIR
r1=/manitou/pmg/projects/korem_lab/Projects/Gherig/mmmbp/tmp/HGF2/${SAMPLE}_1.fastq.gz
r2=/manitou/pmg/projects/korem_lab/Projects/Gherig/mmmbp/tmp/HGF2/${SAMPLE}_2.fastq.gz
CONTIGS=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MMMBPEXT/megahit/${SAMPLE}/final.contigs.fa

bowtie2-build --threads 8 $CONTIGS $TMPDIR/idx #build the index

echo $r1 $r2
bowtie2 --threads 8 -x ${TMPDIR}/idx -U $r1 | samtools view -bS - | samtools sort - -o ${TMPDIR}/alignment_1.bam #align first reads
bowtie2 --threads 8 -x ${TMPDIR}/idx -U $r2 | samtools view -bS - | samtools sort - -o ${TMPDIR}/alignment_2.bam #align second reads
samtools merge -@ 8 ${TMPDIR}/alignment_total.bam ${TMPDIR}/alignment_1.bam ${TMPDIR}/alignment_2.bam #merge the alignments 
samtools sort -@ 8 -n ${TMPDIR}/alignment_total.bam -o ${TMPDIR}/alignment.bam #sort by read names 


/usr/bin/time -v python /burg/pmg/users/ic2465/MetaCarvel/run.py -a $CONTIGS -m ${TMPDIR}/alignment.bam -d $TMPDIR

python /burg/pmg/users/ic2465/Projects/MANU_copangraph/Utils/add_seq_to_mcvl_gfa.py ${TMPDIR}/scaffold_graph.gfa ${CONTIGS}

mkdir -p ${OUTDIR}/${SAMPLE}
mv ${TMPDIR}/*.gfa ${OUTDIR}/${SAMPLE}
mv ${TMPDIR}/*.txt ${OUTDIR}/${SAMPLE}
mv ${TMPDIR}/alignment.bam ${OUTDIR}/${SAMPLE}

