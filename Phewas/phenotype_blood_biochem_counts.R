setwd("/rds/projects/g/gkoutosg-variant-prediction/Yuxi/PheWAS")

library(data.table)

library(readr)

phenotype <- as.data.frame(fread("phenotype.csv"))

phenotype_colname <- colnames(phenotype)

## 1. add covariate phenotaypes!!!

covariates <- read_tsv("covariates.tsv")

covariates_colname <- colnames(covariates)

match_covariates_colomn <- match(covariates_colname,phenotype_colname)

## 2. add blood biochem and counts phenotypes!!!

all_feild_ID <- gsub("x","",sapply(strsplit(phenotype_colname,"_"),'[', 1))

blood_biochem_counts <- read.csv("blood_biochem_counts_feildID.csv")

blood_biochem_counts_feild_ID <- as.character(blood_biochem_counts$Field_ID)

match_blood_biochem_counts_colomn <- match(blood_biochem_counts_feild_ID,all_feild_ID)

# 3. extract selected colomns!!!

phenotype_blood_biochem_counts <- phenotype[,c(match_covariates_colomn,match_blood_biochem_counts_colomn)]

#After generate the file, you need to delete last row "xeid=6026172"

phenotype_blood_biochem_counts1 <- phenotype_blood_biochem_counts[phenotype_blood_biochem_counts$xeid!=6026172,]

save(phenotype_blood_biochem_counts1,file = "phenotype_blood_biochem_counts.RData")

write_tsv(phenotype_blood_biochem_counts1,"phenotype_blood_biochem_counts.tsv")

#Rscript phenotype_blood_biochem_counts.R

