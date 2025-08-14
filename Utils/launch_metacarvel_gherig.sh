#!/bin/bash 
#SBATCH --job-name=mcvlco
#SBATCH --time=120:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/co-mcvl_%x-%j.log
OUTDIR=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/MCVL_${SAMPLE}_CO

TMPDIR=/pmglocal/tk2829/${SAMPLE}_MCVL_CO
rm -rf $TMPDIR
mkdir $TMPDIR

r1=${TMPDIR}/pooled_1.fastq
r2=${TMPDIR}/pooled_2.fastq
cat /burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/${SAMPLE}_sample_names.csv | while read s; do
	zcat /manitou/pmg/projects/korem_lab/Projects/Gherig/mmmbp/tmp/HGF2/${s}_1.fastq.gz >> $r1
	zcat /manitou/pmg/projects/korem_lab/Projects/Gherig/mmmbp/tmp/HGF2/${s}_2.fastq.gz >> $r2
done

CONTIGS=/burg/pmg/users/ic2465/Projects/MANU_copangraph/data/GherigGraphQuality/MAG_results/${SAMPLE}_sample_coasm/final.contigs.fa

bowtie2-build --threads 4 $CONTIGS $TMPDIR/idx #build the index

echo $r1 $r2
bowtie2 --threads 4 -x ${TMPDIR}/idx -U $r1 | samtools view -bS - | samtools sort - -o ${TMPDIR}/alignment_1.bam #align first reads
bowtie2 --threads 4 -x ${TMPDIR}/idx -U $r2 | samtools view -bS - | samtools sort - -o ${TMPDIR}/alignment_2.bam #align second reads
samtools merge -@ 4 ${TMPDIR}/alignment_total.bam ${TMPDIR}/alignment_1.bam ${TMPDIR}/alignment_2.bam #merge the alignments 
samtools sort -@ 4 -n ${TMPDIR}/alignment_total.bam -o ${TMPDIR}/alignment.bam #sort by read names 


/usr/bin/time -v python /burg/pmg/users/ic2465/MetaCarvel/run.py -a $CONTIGS -m ${TMPDIR}/alignment.bam -d $TMPDIR

python /burg/pmg/users/ic2465/Projects/MANU_copangraph/Utils/add_seq_to_mcvl_gfa.py ${TMPDIR}/scaffold_graph.gfa ${CONTIGS}

mkdir -p ${OUTDIR}
mv ${TMPDIR}/*.gfa ${OUTDIR}
mv ${TMPDIR}/*.txt ${OUTDIR}
mv ${TMPDIR}/alignment.bam ${OUTDIR}

