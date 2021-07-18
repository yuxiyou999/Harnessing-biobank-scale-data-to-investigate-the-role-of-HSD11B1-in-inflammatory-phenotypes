setwd("/rds/projects/g/gkoutosg-variant-prediction/Yuxi/PheWAS/results_allSNPs")

# Phewas manhattan plot

library(data.table)
library(qqman)
library(ggpubr)

file_names <- gsub(".csv","",list.files(pattern = "*.csv"))

phewas_outcometable <- list()
PheWasmanPlot <- list()
outcome_info <- as.data.frame(fread("/rds/projects/g/gkoutosg-variant-prediction/Yuxi/PheWAS/PHESANT/variable-info/outcome-info.tsv"))

for (file_number in c(1:length(file_names))) {
  
  file_dir <- paste("/rds/projects/g/gkoutosg-variant-prediction/Yuxi/PheWAS/results_allSNPs/",file_names[file_number],".csv",sep = "")
  
  temp_table <- as.data.frame(fread(file_dir))
  colnames(temp_table)[1] <- "FieldID"
  
  temp_table_withdescription <- left_join(temp_table,outcome_info,by="FieldID")
  temp_table_withdescription$phenotype <- temp_table_withdescription$FieldID
  
  path_split <- strsplit(temp_table_withdescription$Path,">")
  path_colomn <- data.frame(c(1:length(path_split)))
  for (i in c(1:length(path_split))) {
    fullpath <- path_split[[i]]
    path <- fullpath[length(fullpath)]
    path_colomn[i,1] <- path
  }
  temp_table_withdescription$description <- path_colomn[,1]
  temp_table_withdescription<-temp_table_withdescription[temp_table_withdescription$description==" Blood biochemistry" | temp_table_withdescription$description== " Blood count",]
  
  phewas_outcometable[[file_number]] <- temp_table_withdescription

}

# rbind results from all SNPs and extract significant traits
significant_traits <- as.data.frame(t(c(1:dim(phewas_outcometable[[file_number]])[2])))
colnames(significant_traits) <- colnames(phewas_outcometable[[file_number]])
significant_traits$SNP <- ""

#results_from_all_SNPs <- significant_traits

for (file_number in c(1:length(file_names))) {
  
  temp_outcometable <- phewas_outcometable[[file_number]] 
  
  SNP_name <- strsplit(file_names[file_number],"x")[[1]][2]
  
  if (substr(SNP_name,1,2) == "rs") {
    temp_outcometable$SNP <- strsplit(SNP_name,"_")[[1]][1]
  } else {
    temp_outcometable$SNP <- paste(strsplit(SNP_name,"_")[[1]][1],strsplit(SNP_name,"_")[[1]][2],sep="_")
  }
  
  Phewas_full_results <- rbind(significant_traits,temp_outcometable)
  
}



Phewas_full_results <- Phewas_full_results[-1,c(19,17,16,3,20,7)]
colnames(Phewas_full_results)[1] <- "PhenotypeDescription"

write_tsv(Phewas_full_results,"Phewas_full_results.tsv")
#Phewas_full_results <- read_tsv("Phewas_full_results.tsv")

Phewas_full_manplot <- ggplot(Phewas_full_results, aes(y = -log10(pvalue), x = SNP, color = Field)) + 
  geom_point() + 
  geom_hline(yintercept=-log10(5E-8), colour="#990000", linetype="dashed") +
  xlab("3061 SNPs") +
  ylab("-log P value") + 
  theme_bw() +
  theme(legend.margin = margin(0, 0, 0, 0),
        legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(col = "Phenotype\n\nCategory") + 
  guides(colour = guide_legend(nrow = 16,ncol =4)) 

options(bitmapType='cairo')
ggsave(Phewas_full_manplot,filename = "Phewas_full_manplot.png",dpi = 450, width = 12, height = 6)

