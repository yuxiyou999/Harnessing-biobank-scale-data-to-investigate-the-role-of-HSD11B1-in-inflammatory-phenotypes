setwd("/rds/projects/g/gkoutosg-variant-prediction/Yuxi/GWAS/ukb_dataset")
allsampleid <- read.table("allsampleid.txt")

set.seed(123)
index <- sample(2, nrow(allsampleid), replace = T, prob = c(0.5, 0.5))
table(index)
subset1_idlist <- allsampleid[index==1,]
subset2_idlist <- allsampleid[index==2,]

write.table(subset1_idlist,file = "subset1_idlist.txt",row.names = F,col.names = F)
write.table(subset2_idlist,file = "subset2_idlist.txt",row.names = F,col.names = F)


