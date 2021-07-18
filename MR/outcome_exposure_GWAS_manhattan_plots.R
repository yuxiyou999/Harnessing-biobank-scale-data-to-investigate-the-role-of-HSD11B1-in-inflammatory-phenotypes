# file_list,Exposure_list,Outcome_list were come from the MR.R script

Exposure_data_sets <- file_list[file_list %in% Exposure_list$traits]
Outcome_data_sets <- file_list[!(file_list %in% Exposure_list$traits)]

##################

GWAS_Outcome_manPlot <- list()
for (i in c(1:length(Outcome_data_sets))) {
  
  datOutcome <- get(Outcome_data_sets[i])
  
  GWAS_Outcome_title_split <- strsplit(Outcome_data_sets[i],"_")[[1]]
  if (length(GWAS_Outcome_title_split)>5 | nchar(Outcome_data_sets[i])>30) {
    GWAS_Outcome_title <- paste(paste(GWAS_Outcome_title_split[1:3],collapse = " "),"\n",paste(GWAS_Outcome_title_split[4:length(GWAS_Outcome_title_split)],collapse = " "),sep = "")
  } else {
    GWAS_Outcome_title <- paste(gsub("_"," ",Outcome_data_sets[i]),"\n",sep = "")
  }
  
  GWAS_Outcome_manPlot[[i]] <- ggplot(datOutcome, aes(y = LOG10P, x = ID)) +
                                            geom_point() +
                                            geom_hline(yintercept=-log10(5E-8), colour="#990000", linetype="dashed") +
                                            xlab("142 SNPs") +
                                            ylab("-log p") +
                                            ggtitle(GWAS_Outcome_title)  +
                                            theme_bw() +
                                            theme(#panel.border = element_blank(),
                                              axis.text.x=element_blank(),
                                              #axis.ticks.x=element_blank(),
                                              #panel.grid.major.x = element_blank(),
                                              panel.grid.minor.x = element_blank(),
                                              plot.title = element_text(size = 10),
                                              axis.title.x = element_text(size = 10),
                                              axis.title.y = element_text(size = 10)) 

}


GWAS_Outcome_manPlot_all <- ggarrange(GWAS_Outcome_manPlot[[1]],GWAS_Outcome_manPlot[[2]],GWAS_Outcome_manPlot[[3]],
          GWAS_Outcome_manPlot[[4]],GWAS_Outcome_manPlot[[5]],GWAS_Outcome_manPlot[[6]],
          GWAS_Outcome_manPlot[[7]],GWAS_Outcome_manPlot[[8]],GWAS_Outcome_manPlot[[9]],
          ncol = 4, nrow = 3, labels = LETTERS[1:9])
options(bitmapType='cairo')
ggsave(GWAS_Outcome_manPlot_all,filename = "GWAS_Outcome_manPlot_all.png",dpi = 450, width = 12, height = 6)


##################

GWAS_Exposure_manPlot <- list()
for (j in c(1:length(Exposure_data_sets))) {
  
  datExposure <- get(Exposure_data_sets[j])
  
  GWAS_Exposure_title <- gsub("_"," ",Exposure_data_sets[j])
 
  GWAS_Exposure_manPlot[[j]] <- ggplot(datExposure, aes(y = LOG10P, x = ID)) +
    geom_point() +
    geom_hline(yintercept=-log10(5E-8), colour="#990000", linetype="dashed") +
    xlab("142 SNPs") +
    ylab("-log p") +
    ggtitle(GWAS_Exposure_title)  +
    theme_bw() +
    theme(#panel.border = element_blank(),
      axis.text.x=element_blank(),
      #axis.ticks.x=element_blank(),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(size = 10),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)) 
  
}

GWAS_Exposure_manPlot_all <- ggarrange(GWAS_Exposure_manPlot[[1]],GWAS_Exposure_manPlot[[2]],GWAS_Exposure_manPlot[[3]],
                                       GWAS_Exposure_manPlot[[4]],GWAS_Exposure_manPlot[[5]],GWAS_Exposure_manPlot[[6]],
                                       GWAS_Exposure_manPlot[[7]],GWAS_Exposure_manPlot[[8]],GWAS_Exposure_manPlot[[9]],
                                       GWAS_Exposure_manPlot[[10]],GWAS_Exposure_manPlot[[11]],
                                       ncol = 4, nrow = 3, labels = LETTERS[1:length(Exposure_data_sets)])
options(bitmapType='cairo')
ggsave(GWAS_Exposure_manPlot_all,filename = "GWAS_Exposure_manPlot_all.png",dpi = 450, width = 12, height = 6)

