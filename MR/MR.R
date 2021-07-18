setwd("/rds/projects/g/gkoutosg-variant-prediction/Yuxi/MR")
#remotes::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#install.packages("MendelianRandomization")
library(MendelianRandomization)
library(Hmisc) # make the variable's first letter upper case

# step 1: read in data and format
file_name <- list.files(pattern = "*.regenie")

file_list <- rep("XXX",length(file_name))
for (x in c(1:length(file_name))) {
  # prep the name of traits from .regenie file
  file_list[x] <- gsub("phenotypeLINEAR_","",gsub("phenotype_","",strsplit(gsub(".regenie","",file_name[x]),"_chr1_")[[1]][1]))
}
file_list <- capitalize(file_list)


# generate expourse and outcome datasets
for (i in c(1:length(file_name))) {
  
  # generate original expourse or outcome datasets
  trait <- read.table(file_name[i],header = T,sep=" ")
  
  # create pval column
  trait$pval <- 10^(-trait$LOG10P)
  
  # standardize minor allele
  # create dummy variable for saying which allele is minor
  trait$whichMinor <- NA
  # if the allele1 frequency is higher, allele0 is minor
  trait$whichMinor[trait$A1FREQ >= 0.5] <- "ALLELE0"
  # and if lower, it is minor
  trait$whichMinor[trait$A1FREQ < 0.5] <- "ALLELE1"
  table(trait$whichMinor)
  # in this example, we have 8 cases where the minor allele is Allele1, and 92 cases where it is allele0
  tmptraitA1minor <- trait[trait$A1FREQ < 0.5,]
  tmptraitA0minor <- trait[trait$A1FREQ >= 0.5,]
  tmptraitA0minor$ALLELE0n <- tmptraitA0minor$ALLELE1
  tmptraitA0minor$ALLELE1n <- tmptraitA0minor$ALLELE0
  tmptraitA0minor$A1FREQ <- 1 - tmptraitA0minor$A1FREQ
  tmptraitA0minor$ALLELE0 <- tmptraitA0minor$ALLELE0n
  tmptraitA0minor$ALLELE1 <- tmptraitA0minor$ALLELE1n
  tmptraitA0minor$ALLELE0n <- NULL
  tmptraitA0minor$ALLELE1n <- NULL
  trait <- rbind(tmptraitA0minor,tmptraitA1minor)
  rm(tmptraitA0minor,tmptraitA1minor)
  
  assign(file_list[i],trait)
  
}

# the Exposure data set is "Exposure_Mean_corpuscular_volume" etc.
# others are all Outcome data sets!  
print(file_list)

# figure out which are the exposure traits
Exposure_list <- read.csv("Phewas_summary.csv",header = T)
Exposure_list$traits <- gsub(" ","_",Exposure_list$traits)
#Exposure_list$traits %in% file_list

Exposure_data_sets <- file_list[file_list %in% Exposure_list$traits]
Outcome_data_sets <- file_list[!(file_list %in% Exposure_list$traits)]

basicMR_result <- data.frame()
datMRHarm_result <- data.frame()
plieotropy_result <- data.frame()
heterogenity_result <- data.frame()
singleSNP_result <- data.frame()
MR_plots_sig <- list()
k <- 1
# MR analysis
for (j in c(1:length(Exposure_data_sets))) {
  
  for (i in c(1:length(Outcome_data_sets))) {
    
    # filter so snps only affect the exposure, not the outcome, based on GWAS pvalue
    # the pvalue of exposure snps should less than 5E-8 at first!!!
    datExposure <- get(Exposure_data_sets[j])[(get(Exposure_data_sets[j])$pval < 5E-8) & (get(Exposure_data_sets[j])$pval < get(Outcome_data_sets[i])$pval),]
    datOutcome <- get(Outcome_data_sets[i])[get(Outcome_data_sets[i])$ID %in% datExposure$ID,]
    
    # avoid there is no significant snps in exposure data set
    if (dim(datExposure)[1]==0) {
      next
    }
    
    # create a phenotype column, change as needed
    datExposure$Phenotype <- gsub("_"," ",Exposure_data_sets[j])
    datOutcome$Phenotype <- gsub("_"," ",Outcome_data_sets[i])
    
    # format for TwoSampleMR package
    datExposureMR <- format_data(dat = datExposure, type="exposure", phenotype_col = "Phenotype", snp_col = "ID",
                                 beta_col = "BETA", se_col = "SE", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",
                                 other_allele_col = "ALLELE0", pval_col = "pval")
    
    datOutcomeMR <- format_data(dat = datOutcome, type="outcome", phenotype_col = "Phenotype", snp_col = "ID",
                                beta_col = "BETA", se_col = "SE", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1",
                                other_allele_col = "ALLELE0", pval_col = "pval")
    # if a warning message comes up about duplicated SNPs, either ignore or pick which you'd like to use
    
    # create your harmonized data
    datMRHarm <- harmonise_data(
      exposure_dat = datExposureMR, 
      outcome_dat = datOutcomeMR
    )
    
    if (i==1 & j==1) {
      datMRHarm_result <- datMRHarm
    } else {
      datMRHarm_result <- rbind(datMRHarm_result,datMRHarm)
    }
    
    # perform a basic analysis
    # check MR methods
    # command -----> TwoSampleMR::mr_method_list()
    #View(mr_method_list())
    
    # which methods would you like to use? how do you change this?
    # command -----> mr(dat,parameters = default_parameters(),method_list = subset(mr_method_list(), use_by_default)$obj)
    basicMR <- TwoSampleMR::mr(datMRHarm,method_list = subset(mr_method_list(), obj=="mr_ivw")$obj)
   
    if (i==1 & j==1) {
      basicMR_result <- basicMR
    } else {
      basicMR_result <- rbind(basicMR_result,basicMR)
    }
    # View Results; good for a table in supplementary
    # command -----> basicMR
    # which methods were used here? why? what makes them different?
    
    # do a leave-one-out analysis to check for SNPs which dominate the model
    #loo <- TwoSampleMR::mr_singlesnp(datMRHarm)
    # view results
    # command -----> loo
    # so you can make tables out of those
    
    # scatterplot
    if (basicMR$pval*45<0.05) {
      
      ylab_name_split <- strsplit(Outcome_data_sets[i],"_")[[1]]
      if (length(ylab_name_split)>2 & length(ylab_name_split)<5) {
        ylab_name <- paste(paste(ylab_name_split[1:2],collapse = " "),"\n",paste(ylab_name_split[3:length(ylab_name_split)],collapse = " "),sep = "")
      } else if (length(ylab_name_split)>=5) {
        ylab_name <- paste(paste(ylab_name_split[1:2],collapse = " "),"\n",paste(ylab_name_split[3:4],collapse = " "),"\n",paste(ylab_name_split[5:length(ylab_name_split)],collapse = " "),sep = "")
      } else {
        ylab_name <- gsub("_"," ",Outcome_data_sets[i])
      }
      
      MR_plots_sig[[k]] <- TwoSampleMR::mr_scatter_plot(mr_results = basicMR,dat = datMRHarm)[[1]] + xlab(gsub("_"," ",Exposure_data_sets[j])) + 
                        ylab(ylab_name) + theme(axis.title.x = element_text(size = 8.5),axis.title.y = element_text(size = 8.5))
      k <- k+1
      # pdf(paste("MR_plots_sig/","MRScatterPlot_",Exposure_data_sets[j],"_VS_",Outcome_data_sets[i],".pdf",sep = ""))
      #   print(TwoSampleMR::mr_scatter_plot(mr_results = basicMR,dat = datMRHarm))
      # dev.off()
      # 
    }
     
    # forest plot
    p2 <- mr_forest_plot(mr_singlesnp(datMRHarm))

    pdf(paste("MR_plots/","MRForestPlot_",Exposure_data_sets[j],"_VS_",Outcome_data_sets[i],".pdf",sep = ""))
      print(p2[[1]])
    dev.off()

    # now some other analyses:

    # # plieotropy (testing the interept of the MR Egger method)
    plt<-mr_pleiotropy_test(datMRHarm)

    if (i==1 & j==1) {
      plieotropy_result <- plt
    } else {
      plieotropy_result <- rbind(plieotropy_result,plt)
    }
    #
    # heterogenity
    het<-mr_heterogeneity(datMRHarm)

    if (i==1 & j==1) {
      heterogenity_result <- het
    } else {
      heterogenity_result <- rbind(heterogenity_result,het)
    }

    # and looking at just single SNPs
    sin<-mr_singlesnp(datMRHarm)

    if (i==1 & j==1) {
      singleSNP_result <- sin
    } else {
      singleSNP_result <- rbind(singleSNP_result,sin)
    }
    
    
  }

}



all_MR_plots1 <-ggarrange(MR_plots_sig[[12]],MR_plots_sig[[13]],MR_plots_sig[[14]],MR_plots_sig[[15]],
                          ncol = 2, nrow = 2,
                         labels=LETTERS[1:4], font.label = list(size = 12, face = "plain", color ="black"),
                         common.legend = TRUE, legend = "bottom")  


all_MR_plots2 <-ggarrange(MR_plots_sig[[16]],MR_plots_sig[[17]],MR_plots_sig[[18]],MR_plots_sig[[19]],MR_plots_sig[[20]],MR_plots_sig[[21]],
                          ncol = 2, nrow = 3,
                          labels=LETTERS[1:6], font.label = list(size = 12, face = "plain", color ="black"),
                          common.legend = TRUE, legend = "bottom")  


options(bitmapType='cairo')
ggsave(all_MR_plots1,filename = "all_MR_plots1.png",dpi = 450, width = 12, height = 6)
ggsave(all_MR_plots2,filename = "all_MR_plots2.png",dpi = 450, width = 12, height = 6)


# bonferroni method adjusts P value!

basicMR_result$bonferroni_pval <- p.adjust(basicMR_result$pval, method = "bonferroni")

write.csv(basicMR_result,"basicMR_result.csv")
write.csv(datMRHarm_result,"datMRHarm_result.csv")
write.csv(plieotropy_result,"plieotropy_result.csv")
write.csv(heterogenity_result,"heterogenity_result.csv")
write.csv(singleSNP_result,"singleSNP_result.csv")

