#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --qos castles
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 40:0:0
#SBATCH --mem 25G
#SBATCH --job-name binary_step2_phenotype_crohn_disease_regional_enteritis
module purge; module load bluebear

#inputs are the phenoFile, covarFile and pred file generated from previous step. Results are output into results directory.

/rds/projects/k/karwatha-rate-af/datacentre/GWAS/regenie/regenie \
    --step 2 \
    --bgen /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/chr1_HSD11B1.bgen --with-bgi \
    --sample /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/chr1_HSD11B1.sample \
    --phenoFile /rds/projects/g/gkoutosg-variant-prediction/Yuxi/phenotypic_data/phenotype_crohn_disease_regional_enteritis.phen \
    --covarFile /rds/projects/g/gkoutosg-variant-prediction/Yuxi/phenotypic_data/full_inclusive.cov \
    --keep /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/subset1_idlist.txt \
    --extract snplist_extracted_from_phewas.snplist \
    --bt \
    --firth 0.01 --approx --firth-se \
    --pred phenotype_crohn_disease_regional_enteritis_step1_pred.list \
    --lowmem /rds/projects/k/karwatha-rate-af/datacentre/GWAS/regenie/tmp \
    --bsize 400 \
    --split \
    --threads 5 \
    --out results/phenotype_crohn_disease_regional_enteritis_chr1_HSD11B1_step2
