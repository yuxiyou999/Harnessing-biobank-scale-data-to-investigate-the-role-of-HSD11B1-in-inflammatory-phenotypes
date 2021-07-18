getwd()
setwd("/rds/projects/g/gkoutosg-variant-prediction/Yuxi/phenotypic_data")

library(readr)
library(dplyr)
library(tidyverse)

############################
### rheumatoid arthritis ###
############################

## grep seropositive and other rheumatoid arthritis 

rheumatoid_arthritis_func <- function(x,pos) {
  x <- x[,c("f.eid","f.131848.0.0","f.131850.0.0")]
  # "f.131848.0.0" ---> "seropositive rheumatoid arthritis"
  # "f.131850.0.0" ---> "other rheumatoid arthritis"
  return(x)
}

rheumatoid_arthritis <- readr::read_delim_chunked("ukb44682.diagnosisdates.tab",delim = "\t",col_names = TRUE,DataFrameCallback$new(rheumatoid_arthritis_func))

table(rheumatoid_arthritis$f.131850.0.0)
summary(rheumatoid_arthritis)

## determine if rheumatoid arthritis diagnosis = 0 or 1

rheumatoid_arthritis$FID <- rheumatoid_arthritis$f.eid
rheumatoid_arthritis$IID <- rheumatoid_arthritis$f.eid

index <- is.na(rheumatoid_arthritis$f.131850.0.0) & is.na(rheumatoid_arthritis$f.131848.0.0)
#table(index)
rheumatoid_arthritis$rheumatoid_arthritis_diagnosis <- NA
rheumatoid_arthritis$rheumatoid_arthritis_diagnosis[index == TRUE] <- 0
rheumatoid_arthritis$rheumatoid_arthritis_diagnosis[index == FALSE] <- 1
#hist(rheumatoid_arthritis$rheumatoid_arthritis_diagnosis)
#table(rheumatoid_arthritis$rheumatoid_arthritis_diagnosis)

phenotype_rheumatoid_arthritis <- rheumatoid_arthritis[,c(4:6)]
write.table(phenotype_rheumatoid_arthritis, file = "phenotype_rheumatoid_arthritis.phen", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

############################
###  C-reactive protein  ###
############################

## grep C-reactive protein colomn 

CRP_func <- function(x,pos) {
  x <- x[,c("f.eid","f.30710.0.0","f.30710.1.0")]
  #field.showcase	field.html	field.tab	col.type	col.name
  #eid	eid	f.eid	Sequence	eid
  #...
  #30710	30710-0.0	f.30710.0.0	Continuous	creactive_protein_f30710_0_0 ---> instance 0
  #30710	30710-1.0	f.30710.1.0	Continuous	creactive_protein_f30710_1_0 ---> instance 1
  #...
  return(x)
}

CRP <- readr::read_delim_chunked("ukb46518.tab",delim = "\t",col_names = TRUE,DataFrameCallback$new(CRP_func))

CRP$FID <- CRP$f.eid
CRP$IID <- CRP$f.eid
CRP$CRP_value <- CRP$f.30710.0.0
CRP$log_CRP_value <- log(CRP$f.30710.0.0)
#hist(CRP$CRP_value)
#hist(CRP$log_CRP_value)

phenotype_CRP <- CRP[,c(4:6)]
write.table(phenotype_CRP, file = "phenotypeLINEAR_CRP.phen", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

### the following codes are written on 2nd July.

##################################
###  extract other phenotypes  ###   
##################################

# Data-Field 	disease
# 131492	Date J44 first reported (other chronic obstructive pulmonary disease)
# 131894	Date M32 first reported (systemic lupus erythematosus)
# 131854	Date M08 first reported (juvenile arthritis)
# 131912	Date M45 first reported (ankylosing spondylitis)
# 131626	Date K50 first reported (crohn's disease [regional enteritis])
# 131852	Date M07 first reported (psoriatic and enteropathic arthropathies)
# 131628	Date K51 first reported (ulcerative colitis)
# 131832	Date L95 first reported (vasculitis limited to skin, not elsewhere classified)

trait_list <- c("other_chronic_obstructive_pulmonary_disease","systemic_lupus_erythematosus",
                "juvenile_arthritis","ankylosing_spondylitis","crohn_disease_regional_enteritis",
                "psoriatic_and_enteropathic_arthropathies","ulcerative_colitis","vasculitis_limited_to_skin")

feid_list <- c("f.131492.0.0","f.131894.0.0","f.131854.0.0","f.131912.0.0","f.131626.0.0","f.131852.0.0","f.131628.0.0","f.131832.0.0")

extract_traits_func <- function(x,pos) {
  x <- x[,c("f.eid",feid_list)]
  return(x)
}

extract_traits <- readr::read_delim_chunked("ukb44682.diagnosisdates.tab",delim = "\t",col_names = TRUE,DataFrameCallback$new(extract_traits_func))

# generate .phen file

for (i in c(1:length(trait_list))) {
  
  ## determine if diagnosis = 0 or 1
  
  trait <- extract_traits[,c(1,(i+1))]
  trait$FID <- trait$f.eid
  trait$IID <- trait$f.eid
  
  index <- is.na(trait[,2]) 
  #table(index)
  trait$trait_diagnosis <- NA
  trait$trait_diagnosis[index == TRUE] <- 0
  trait$trait_diagnosis[index == FALSE] <- 1
  #table(trait$trait_diagnosis)
  
  phenotype <- trait[,c(3:5)]
  write.table(phenotype, file = paste("phenotype_",trait_list[i],".phen",sep = ""), sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}

# "juvenile_arthritis" contains no 1, which means there is no diagnosed case!


#################################
###  Mean corpuscular volume  ###
#################################
# Data-Field 30040
# Description:	Mean corpuscular volume
# Category:	Blood count - Blood assays - Biological samples
# Instance: 0, 1, 2

## grep Mean corpuscular volume colomn

Mean_corpuscular_volume_func <- function(x,pos) {
  x <- x[,c("f.eid","f.30040.0.0")]
  return(x)
}

Mean_corpuscular_volume <- readr::read_delim_chunked("ukb46518.tab",delim = "\t",col_names = TRUE,DataFrameCallback$new(Mean_corpuscular_volume_func))

Mean_corpuscular_volume$FID <- Mean_corpuscular_volume$f.eid
Mean_corpuscular_volume$IID <- Mean_corpuscular_volume$f.eid
Mean_corpuscular_volume$Mean_corpuscular_volume_value <- Mean_corpuscular_volume$f.30710.0.0

#hist(Mean_corpuscular_volume$Mean_corpuscular_volume_value)

phenotype_Mean_corpuscular_volume <- Mean_corpuscular_volume[,c(3:5)]
write.table(phenotype_Mean_corpuscular_volume, file = "phenotypeLINEAR_Mean_corpuscular_volume.phen", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)



########################################
###  extract other blood phenotypes  ###   
########################################

# Data-Field	traits	number of snps
# 30050	Mean corpuscular haemoglobin	111
# 30040	Mean corpuscular volume	106
# 30260	Mean reticulocyte volume	106
# 30010	Red blood cell (erythrocyte) count	106
# 30770	IGF-1	44
# 30270	Mean sphered cell volume	37
# 30720	Cystatin C	19
# 30090	Platelet crit	12
# 30880	Urate	8
# 30140	Neutrophill count	6
# 30650	Aspartate aminotransferase	4

Phewas_summary <- read.csv("Phewas_summary.csv")
# delete the .phen file has generated : "Mean corpuscular volume"
Phewas_summary <- Phewas_summary[Phewas_summary$traits!="Mean corpuscular volume",]
Phewas_summary$feid <- paste("f.",Phewas_summary$Data.Field,".0.0",sep = "")

trait_list <- Phewas_summary$traits
feid_list <- Phewas_summary$feid

extract_traits_func <- function(x,pos) {
  x <- x[,c("f.eid",feid_list)]
  return(x)
}

extract_traits <- readr::read_delim_chunked("ukb46518.tab",delim = "\t",col_names = TRUE,DataFrameCallback$new(extract_traits_func))
#write_csv(extract_traits,"extract_traits.csv",col_names = T)
#extract_traits <- read_csv("extract_traits.csv",col_names = T)

# generate .phen file

for (i in c(1:length(trait_list))) {
  
  trait <- as.data.frame(extract_traits[,c(1,(i+1))])
  trait$FID <- trait$f.eid
  trait$IID <- trait$f.eid
  trait$trait_value <- trait[,2]
  
  phenotype <- trait[,c(3:5)]
  write.table(phenotype, file = paste("phenotypeLINEAR_",gsub(" ","_",trait_list[i]),".phen",sep = ""), sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}
