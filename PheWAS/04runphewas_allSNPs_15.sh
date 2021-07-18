#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --qos castles
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --array 1-100
#SBATCH --time 5:0:0
#SBATCH --mem 20G
#SBATCH --job-name trial-run
module purge; module load bluebear
module load R/3.6.0-foss-2019a

#this runs the phewas, across each snp per job
#the phewas is run within the PHESANT folder so everything is routed back to the previous working directory

snp=$(head -n1 alt/extract_chr1_13_allSNPs_15.csv | tr ',' '\n' | tail -n+2 | sed "${SLURM_ARRAY_TASK_ID}q;d")

pt=$1

mkdir -p results_allSNPs/${snp}
wd=$(pwd)

cd PHESANT/WAS 

Rscript phenomeScan.r \
--phenofile="${wd}/phenotype_blood_biochem_counts.tsv" \
--traitofinterestfile="${wd}/alt/extract_chr1_13_allSNPs_15.csv" \
--variablelistfile="../variable-info/outcome-info.tsv" \
--datacodingfile="../variable-info/data-coding-ordinal-info.txt" \
--traitofinterest="${snp}" \
--resDir="${wd}/results_allSNPs/${snp}/" \
--sensitivity \
--genetic=TRUE \
--userId="xeid" \
--tab=TRUE 