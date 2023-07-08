# GSE192714  Mus
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9630274/
# Result -- Mll3 or Mll4 ablation induces transcriptional priming of the GSDMD-mediated pyroptotic pathway
# indicating that Mll3 and Mll4 deletion induces GSDMD expression and could transcriptionally prime tumor cells for pyroptosis. 


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
library(dplyr)


# ==============================================================================
# GEO Datasets prepare
# ==============================================================================

# GSE192714_pdata
GSE192714_pData <- read.csv("../../Data/PScore_validation/GSE192714_series_matrix.txt.gz",header = FALSE,sep = "\t",col.names = paste0("V",1:15))
GSE192714_pData <- GSE192714_pData[substr(GSE192714_pData$V1,1,8)=="!Sample_",]
GSE192714_pData$V1 <- make.names(gsub("!Sample_","",GSE192714_pData$V1),unique=T)
row.names(GSE192714_pData) <- GSE192714_pData[,1]
GSE192714_pData <- GSE192714_pData[,-1]
GSE192714_pData <- as.data.frame(t(GSE192714_pData))
row.names(GSE192714_pData) <- GSE192714_pData$geo_accession
GSE192714_pData <- GSE192714_pData[1:8,]

# GSE192714_matrix
options(scipen=200)
expr_txt <- list.files("../../Data/PScore_validation/GSE192714_RAW/")
GSE192714_count <-  lapply(expr_txt,function(x)read.table(paste0("../../Data/PScore_validation/GSE192714_RAW/",x),sep="\t",header = T))
names(GSE192714_count) <- expr_txt

# merge samples
tmp_1 <- GSE192714_count[[1]]
tmp_1$Ensembl_Gene_ID <- gsub("\\..*$","",tmp_1$Ensembl_Gene_ID)
for(i in 2:length(GSE192714_count)){
  tmp_2 <- GSE192714_count[[i]]
  tmp_2$Ensembl_Gene_ID <- gsub("\\..*$","",tmp_2$Ensembl_Gene_ID)
  tmp_1 <- merge(tmp_1,tmp_2,by = "Ensembl_Gene_ID")
}
GSE192714_count_merge <- tmp_1

# exon length (count to fpkm)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
exons_gene <- exonsBy(txdb, by = "gene") 
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_lens[1:10]
gene_length <- sapply(exons_gene_lens,function(x){x})
id_length <- as.data.frame(gene_length)
id_length$Ensembl_Gene_ID <- rownames(id_length)
save(id_length,file = "d:/Work/Data/genome/mouse/GRCm38(mm10)/gene_length_TxDb.Mmusculus.UCSC.mm10.ensGene.RData")

GSE192714_count_len <- merge(id_length, GSE192714_count_merge, by = "Ensembl_Gene_ID")
  
# id convert
library(org.Mm.eg.db)
gene.tmp <- select(org.Mm.eg.db, keys=GSE192714_count_len$Ensembl_Gene_ID, columns=c("SYMBOL"), keytype="ENSEMBL" )
gene.tmp <- gene.tmp[!duplicated(gene.tmp$ENSEMBL),]
GSE192714_symbol <- merge(gene.tmp, GSE192714_count_len, by.x = "ENSEMBL", by.y = "Ensembl_Gene_ID")
GSE192714_symbol <- GSE192714_symbol[!duplicated(GSE192714_symbol$SYMBOL),]
GSE192714_symbol <- GSE192714_symbol[!is.na(GSE192714_symbol$SYMBOL),]
# [1] 27567    11

# countToFpkm  mm10(GSE192714_pData$data_processing.3)
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )}

GSE192714_FPKM <- as.data.frame(apply(GSE192714_symbol[4:11],2,countToFpkm,effLen = GSE192714_symbol$gene_length))
rownames(GSE192714_FPKM) <- GSE192714_symbol$SYMBOL
GSE192714_FPKM <- as.matrix(GSE192714_FPKM)

# > dim(GSE192714_FPKM)
# [1] 26560     8
# GSE192714_FPKM[1:4,1:4]
# Count_WT_Rep1 Count_WT_Rep2 Count_WT_Rep3 Count_WT_Rep4
# Gnai3   18.96909645    21.1308760      16.06056      17.13688
# Pbsn     0.00000000     0.0000000       0.00000       0.00000
# H19      0.00000000     0.1231331       0.00000       0.00000
# Scml2    0.04917316     0.0000000       0.00000       0.00000


# ==============================================================================
# score
# ==============================================================================
library(homologene)

load("coxForest/uniCoxSig_list.RData")
interest_genes <- uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"]

interest_genes_mus <- homologene(interest_genes, inTax = 9606, outTax = 10090)
interest_genes[interest_genes %in% interest_genes_mus$`9606` == FALSE]
# [1] "CASP5"  "CHMP4A" "GSDMB"  "GZMB"   "NLRP7" 

genelist <- list(PS=c(interest_genes_mus$`10090`,"Gzmb"))
genelist$PS[genelist$PS %in% rownames(GSE192714_FPKM) ==FALSE]
# "Dfnb59" 
c("Dfnb59", "Gm1001", "Pjvk") %in% rownames(GSE192714_FPKM)
rownames(GSE192714_FPKM)[rownames(GSE192714_FPKM)=="Pjvk"] <- "Dfnb59"

Pyro_matrix <- as.matrix(GSE192714_FPKM)
data_ES <- GSVA::gsva(Pyro_matrix,genelist,method="ssgsea")
data_ES <- t(data_ES) %>% data.frame()
head(data_ES)
data_ES_zscore <- scale(data_ES) %>% data.frame()



# ==============================================================================
# visualize
# ==============================================================================
library(ggplot2)
library(ggpubr) 

data_ES_zscore$Group <- rownames(data_ES_zscore)
data_ES_zscore$group <- gsub("Count_|Count_Mll[3-4]{1}|_Rep[1-4]{1}","",data_ES_zscore$Group)
data_ES_zscore$group <- factor(data_ES_zscore$group,levels = c("WT","KO"))
table(data_ES_zscore$group)


p1 <- ggplot(data_ES_zscore,aes(x=group,y=PS,color=group))+
  geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+
  ggbeeswarm::geom_beeswarm(alpha=0.5,size=0.5)+
  stat_compare_means(method = "wilcox.test") +
  theme_bw()
p1


# save data
GSE192714_ggplot <- data_ES_zscore[,c("group","Group","PS")]
GSE192714_ggplot$geo_accession <- "GSE192714"
GSE192714_ggplot$note <- "Mus. B16-F10 cells. Lack CASP5,CHMP4A,GSDMB,NLRP7 which without homologous gene in hsa"
save(GSE192714_ggplot,file="PScoreValidation/GSE192714_ggplot.RData")


