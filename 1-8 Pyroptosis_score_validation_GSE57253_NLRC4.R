# GSE57253


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
library(dplyr)
library(ggplot2)
library(ggpubr) 



# ==============================================================================
# GEO Datasets prepare
# ==============================================================================

# GSE57253_pdata
GSE57253_pData <- read.csv("../../Data/PScore_validation/GSE57253_series_matrix.txt.gz",header = FALSE,sep = "\t",col.names = paste0("V",1:27))
GSE57253_pData <- GSE57253_pData[substr(GSE57253_pData$V1,1,8)=="!Sample_",]
GSE57253_pData$V1 <- make.names(gsub("!Sample_","",GSE57253_pData$V1),unique=T) 

row.names(GSE57253_pData) <- GSE57253_pData[,1]
GSE57253_pData <- GSE57253_pData[,-1]

GSE57253_pData <- as.data.frame(t(GSE57253_pData))
row.names(GSE57253_pData) <- GSE57253_pData$geo_accession


# GSE57253_matrix
options(scipen=200)
GSE57253_RPKM <- read.csv("../../Data/PScore_validation/GSE57253_140225_140305_130716_RNA-Seq_results.transcript.rpkm.txt.gz",
                          skip=1,header = T,sep = "\t")

GSE57253_RPKM <- GSE57253_RPKM[GSE57253_RPKM$Symbol!="#N/A",] 
GSE57253_RPKM$Symbol <- gsub(" .*$","",GSE57253_RPKM$Symbol)
GSE57253_RPKM <- GSE57253_RPKM[!duplicated(GSE57253_RPKM$Symbol),] 

GSE57253_RPKM[,2:27] <- apply(GSE57253_RPKM[,2:27],2,as.numeric)  
rownames(GSE57253_RPKM) <- GSE57253_RPKM$Symbol

GSE57253_RPKM <- GSE57253_RPKM[,substr(colnames(GSE57253_RPKM),1,1)=="s"] 
GSE57253_RPKM <- as.matrix(GSE57253_RPKM)
GSE57253_RPKM[1:4,1:4]

  

# ==============================================================================
# score
# ==============================================================================
load("coxForest/uniCoxSig_list.RData")

interest_genes <- uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"]
genelist <- list(PS=interest_genes);rm(uniCoxSig_list)

interest_genes[interest_genes %in% rownames(GSE57253_RPKM) ==FALSE]

Pyro_matrix <- GSE57253_RPKM
data_ES <- GSVA::gsva(Pyro_matrix,genelist,method="ssgsea")
data_ES <- t(data_ES) %>% data.frame()
head(data_ES)
data_ES_zscore <- scale(data_ES) %>% data.frame()



# ==============================================================================
# merge data
# ==============================================================================

data_ES_zscore$id <- toupper(rownames(data_ES_zscore))
GSE57253_pData$id <- toupper(GSE57253_pData$title)
data_ES_zscore_merge <- merge(GSE57253_pData,data_ES_zscore,by="id")

table(data_ES_zscore_merge$source_name_ch1)
data_ES_zscore_merge$source_name_ch1[data_ES_zscore_merge$source_name_ch1
                                     =="healthy pediatric controls"] <- "Healthy Controls"
data_ES_zscore_merge$source_name_ch1[data_ES_zscore_merge$source_name_ch1
                                     =="NOMID patients with active disease prior to anakinra treatment"] <- "NOMID"
data_ES_zscore_merge$source_name_ch1[data_ES_zscore_merge$source_name_ch1
                                     =="NOMID patients with inactive disease after anakinra treatment"] <- "NOMID after treatment"
data_ES_zscore_merge$source_name_ch1[data_ES_zscore_merge$source_name_ch1
                                     =="patient with NLRC4-MAS"] <- "NLRC4-MAS"
data_ES_zscore_merge <- data_ES_zscore_merge[data_ES_zscore_merge$source_name_ch1!="NOMID after treatment",]
data_ES_zscore_merge$Group <- factor(data_ES_zscore_merge$source_name_ch1,
                                     levels = c("Healthy Controls","NOMID","NLRC4-MAS")) 

# save data
GSE57253_ggplot <- data_ES_zscore_merge[,c("Group","PS")]
GSE57253_ggplot$group <- GSE57253_ggplot$Group
GSE57253_ggplot$geo_accession <- "GSE57253"
GSE57253_ggplot$note <- "homo sapiens"
save(GSE57253_ggplot,file="PScoreValidation/GSE57253_ggplot.RData")


# ==============================================================================
# visualize 
# ==============================================================================
# pdf("PScoreValidation/GSE57253_rpkm.pdf",width = 7,height = 7)
ggplot(data_ES_zscore_merge,aes(x=Group,y=PS,color=Group))+
  geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+
  ggbeeswarm::geom_beeswarm(alpha=0.5,size=0.5)+
  stat_compare_means(method = "kruskal.test") +
  theme_bw()
# dev.off()




