#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --qos castles
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 5:0:0
#SBATCH --mem 100G
#SBATCH --job-name making_subsets
module purge; module load bluebear
module load PLINK/1.9b_6.17-x86_64

plink --bfile biobank --keep subset1_idlist.txt --make-bed --out subset1
plink --bfile biobank --keep subset2_idlist.txt --make-bed --out subset2
