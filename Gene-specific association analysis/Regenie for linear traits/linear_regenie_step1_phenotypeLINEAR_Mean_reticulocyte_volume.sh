#!/bin/bash
#SBATCH --mail-type NONE
#SBATCH --qos castles
#SBATCH --nodes 1
#SBATCH --ntasks 1 
#SBATCH --cpus-per-task 15 
#SBATCH --time 40:0:0 
#SBATCH --mem 40G 
#SBATCH --job-name linear_step1_phenotypeLINEAR_Mean_reticulocyte_volume
module purge; module load bluebear 
 
#input is the phenoFile and potentially the covarFile. Output is a pred file which tells the next step where to look for a loco file which is the NULL model, name is set --out 
/rds/projects/k/karwatha-rate-af/datacentre/GWAS/regenie/regenie \
	--step 1 \
	--bed /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/subset2 \
	--extract /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/biobank_qcd.snplist \
	--keep /rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset/biobank_qcd.id \
	--phenoFile /rds/projects/g/gkoutosg-variant-prediction/Yuxi/phenotypic_data/phenotypeLINEAR_Mean_reticulocyte_volume.phen \
	--covarFile /rds/projects/g/gkoutosg-variant-prediction/Yuxi/phenotypic_data/full_inclusive.cov \
	--bsize 1000 \
	--lowmem /rds/projects/k/karwatha-rate-af/datacentre/GWAS/regenie/tmp \
	--threads 15 \
	--out phenotypeLINEAR_Mean_reticulocyte_volume_step1

#./run_lin.sh
