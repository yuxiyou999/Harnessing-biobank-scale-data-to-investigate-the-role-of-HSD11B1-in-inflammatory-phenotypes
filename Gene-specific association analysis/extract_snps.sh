#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --qos castles
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 10:0:0
#SBATCH --mem 30G
#SBATCH --job-name trial-run
module purge; module load bluebear
module load QCTOOL/2.0.8-foss-2019b

#this script goes through all of the genetic files and checks if the rsids exist in those files and extracts then
file="ukb_imp_chr1_v3.bgen"
samp="ukb31224_imp_chr1_v3_s487324.sample"

#extract rsids from larger genetic files
qctool -g $file -s $samp \
    -incl-range 01:208859524-210908274 \
    -og chr1_HSD11B1.bgen -os chr1_HSD11B1.sample

#create index
./bgenix -g chr1_HSD11B1.bgen -index

