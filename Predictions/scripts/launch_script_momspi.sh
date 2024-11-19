#!/bin/bash
#SBATCH --job-name=sgc_momspi
#SBATCH --time=48:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --nodelist=m009
#SBATCH --account pmg
#SBATCH --output=/burg/pmg/users/ic2465/Projects/MANU_copangraph/%x-%j.log
wkdir=/pmglocal/ic2465
SGC_DIR=/burg/pmg/users/ic2465/Projects/MANU_copangraph/spacegraphcats


### do stuff

# set up spacegraphcats dir
#mkdir $wkdir
#cd $wkdir
#cp -r $SGC_DIR ./
#cd spacegraphcats/inputs
#mv mgx_raw_momspi mgx_raw
#mv metadata.csv.momspi metadata.csv
#cd ..
#tar -xvf snkdmp.tar
#
#snakemake -s 00_select_query_species_for_dda.snakefile --use-conda --rerun-incomplete -c 32
cd /pmglocal/ic2465/spacegraphcats
snakemake -s 01_perform_dda.snakefile --use-conda --rerun-incomplete -c 16 --unlock
snakemake -s 01_perform_dda.snakefile --use-conda --rerun-incomplete -c 16
