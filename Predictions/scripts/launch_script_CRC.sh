#!/bin/bash
#SBATCH --job-name=sgc_crc
#SBATCH --time=72:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --nodelist=m012
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/%x-%j.log
wkdir=/pmglocal/ic2465
SGC_DIR=/burg/pmg/users/ic2465/Projects/MANU_copangraph/spacegraphcats


### do stuff

# set up spacegraphcats dir
#mkdir $wkdir
#cd $wkdir
#cp -r $SGC_DIR ./spacegraphcats_crc
#cd spacegraphcats_crc/inputs
#for f in /burg/pmg/users/az2732/CRC_CN_HGF2_gzipped/*.fastq.gz; do
#	echo $f
#	cp $f ./mgx_raw_crc
#done
#
#cd mgx_raw_crc
## rename for sgc
#for f in *.fastq.gz; do
#	mv ${f} ${f%.fastq.gz}.fq.gz
#done
## reformat for sgc
#for f in *.fq.gz; do
#	gunzip $f
#	flat_f=${f%.gz}
#	sed -i 's/[[:space:]]/_/g' $flat_f
#	sed -i 's/_length=[0-9]\+$//' $flat_f
#	pigz $flat_f
#done
#cd ..
#mv mgx_raw_crc mgx_raw
#cd ..

# untar snakemake stuff
#tar -xvf snkdmp.tar
cd $wkdir
cd spacegraphcats_crc
#snakemake -s 00_select_query_species_for_dda.snakefile --use-conda --rerun-incomplete -c 32
snakemake -s 01_perform_dda.snakefile --use-conda --rerun-incomplete -c 16 --unlock
snakemake -s 01_perform_dda.snakefile --use-conda --rerun-incomplete -c 16  
