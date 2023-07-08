# GSE153494   	Mus musculus


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
library(dplyr)


# ==============================================================================
# data prepare
# ==============================================================================

# GSE153494_pdata
GSE153494_pData <- read.csv("../../Data/PScore_validation/GSE153494_series_matrix.txt.gz",header = FALSE,sep = "\t",col.names = paste0("V",1:19))
GSE153494_pData <- GSE153494_pData[substr(GSE153494_pData$V1,1,8)=="!Sample_",]
GSE153494_pData$V1 <- make.names(gsub("!Sample_","",GSE153494_pData$V1),unique=T) 

row.names(GSE153494_pData) <- GSE153494_pData[,1]
GSE153494_pData <- GSE153494_pData[,-1]

GSE153494_pData <- as.data.frame(t(GSE153494_pData))
row.names(GSE153494_pData) <- GSE153494_pData$geo_accession


# GSE153494_matrix
options(scipen=200)
GSE153494_RPKM <- read.csv("../../Data/PScore_validation/GSE153494_All_sample_RPKM.csv.gz",header = T,sep = ",")
GSE153494_RPKM <- GSE153494_RPKM[!duplicated(GSE153494_RPKM$X),]
rownames(GSE153494_RPKM) <- GSE153494_RPKM$X

GSE153494_RPKM[,2:19] <- apply(GSE153494_RPKM[,2:19],2,as.numeric)
GSE153494_RPKM <- GSE153494_RPKM[,grep(pattern="FPKM",colnames(GSE153494_RPKM))]

GSE153494_RPKM <- as.matrix(GSE153494_RPKM)
GSE153494_RPKM[1:4,1:4]



# ==============================================================================
# score
# ==============================================================================
load("coxForest/uniCoxSig_list.RData")

interest_genes <- uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"]

library(homologene)
homologene::taxData
# 10090   Mus musculus
# 9606    Homo sapiens

interest_genes_mus <- homologene(interest_genes, inTax = 9606, outTax = 10090)
human2mouse(interest_genes)
genelist <- list(PS=interest_genes_mus$`10090`)


# [1] "Aim2"     "Apip"     "Casp1"    "Casp3"    "Casp4"    "Casp8"    "Cflar"    "Chmp2b"   "Chmp5"    "Dfnb59"   "Gsdmd"   
# [12] "Gzma"     "Il18"     "Il1b"     "Irf1"     "Irf2"     "Mefv"     "Naip2"    "Nlrc4"    "Nlrp1b"   "Nlrp3"    "Nlrp6"   
# [23] "Nod2"     "Tlr4"     "Tnfrsf1b" "Zbp1"    

interest_genes[interest_genes%in%interest_genes_mus$`9606`==FALSE]
# [1] "CASP5"  "CHMP4A" "GSDMB"  "GZMB"   "NLRP7" 
# genelist$PS <- c(genelist$PS,"Gzmb")
tmp <- human2mouse(c("AIM2", "GSDMC", "GSDMD", "IL18", "NLRP6", "PRKACA"))
genelist <- list(#PRGs_74 = PRGs_74$Symbol.ID,
  PRGs_prot =  c(genelist$PS,"Gzmb"),
  Wang_Ann.Transl.Med_2022 = tmp$mouseGene)
genelist$PRGs_prot[genelist$PRGs_prot %in% rownames(GSE153494_RPKM) ==FALSE]


Pyro_matrix <- GSE153494_RPKM
data_ES <- GSVA::gsva(Pyro_matrix,genelist,method="ssgsea")
data_ES <- t(data_ES) %>% data.frame()
head(data_ES)
data_ES_zscore <- scale(data_ES) %>% data.frame()



# ==============================================================================
# merge data
# ==============================================================================

data_ES_zscore$id <- gsub("FPKM.","",rownames(data_ES_zscore))
data_ES_zscore$id <- gsub("Mlh","M1h",data_ES_zscore$id)
data_ES_zscore$id <- gsub("M6n","M6h",data_ES_zscore$id)
GSE153494_pData$id <- gsub(" rep","",GSE153494_pData$title)

data_ES_zscore_merge <- merge(GSE153494_pData,data_ES_zscore,by="id")

data_ES_zscore_merge$group <- gsub("M| rep.*$","",data_ES_zscore_merge$title)
data_ES_zscore_merge$group[data_ES_zscore_merge$group=="con"] <- "control"

table(data_ES_zscore_merge$group)

data_ES_zscore_merge$Group <- factor(data_ES_zscore_merge$group,
                                     levels = c("control","10m","1h","6h","24h","72h")) 

# save data
GSE153494_ggplot <- data_ES_zscore_merge[,c("Group","PS")]
GSE153494_ggplot$group <- factor(GSE153494_ggplot$Group,levels=c("control", "10m", "1h", "6h", "24h", "72h"))
GSE153494_ggplot$geo_accession <- "GSE153494"
GSE153494_ggplot$note <- "mus musculus, lack CASP5, CHMP4A, GSDMB, NLRP7"
save(GSE153494_ggplot,file="PScoreValidation/GSE153494_ggplot.RData")



# ==============================================================================
# visualize
# ==============================================================================
library(ggplot2)
library(ggpubr) 

pdf("PScoreValidation/GSE153494_RPKM_mus.pdf",width = 5,height = 5)
ggplot(data_ES_zscore_merge,aes(x=Group,y=PS,color=Group))+
  geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+
  ggbeeswarm::geom_beeswarm(alpha=0.5,size=0.5)+
  stat_compare_means(method = "kruskal.test") +
  theme_bw()
dev.off()

ggplot(data_ES_zscore_merge,aes(x=Group,y=Wang_Ann.Transl.Med_2022,color=Group))+
  geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+
  ggbeeswarm::geom_beeswarm(alpha=0.5,size=0.5)+
  stat_compare_means(method = "kruskal.test") +
  theme_bw()


