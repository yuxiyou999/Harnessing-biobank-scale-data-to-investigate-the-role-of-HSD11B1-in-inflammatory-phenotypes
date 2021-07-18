#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --qos castles
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5
#SBATCH --time 40:0:0
#SBATCH --mem 25G
#SBATCH --job-name phenotypeLINEAR_Mean_corpuscular_volume_step2
module purge; module load bluebear

#inputs are the phenoFile, covarFile and pred file generated from previous step. Results are output into results directory.

/rds/projects/k/karwatha-rate-af/datacentre/GWAS/regenie/regenie \
    --step 2 \
    --bgen /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/chr1_HSD11B1.bgen --with-bgi \
    --sample /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/chr1_HSD11B1.sample \
    --phenoFile /rds/projects/g/gkoutosg-variant-prediction/Yuxi/phenotypic_data/phenotypeLINEAR_Mean_corpuscular_volume.phen \
    --covarFile /rds/projects/g/gkoutosg-variant-prediction/Yuxi/phenotypic_data/full_inclusive.cov \
    --keep /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/subset2_idlist.txt \
    --extract snplist_extracted_from_phewas.snplist \
    --pred phenotypeLINEAR_Mean_corpuscular_volume_step1_pred.list \
    --lowmem /rds/projects/k/karwatha-rate-af/datacentre/GWAS/regenie/tmp \
    --bsize 400 \
    --split \
    --threads 5 \
    --out results/phenotype_Mean_corpuscular_volume_chr1_HSD11B1_step2


