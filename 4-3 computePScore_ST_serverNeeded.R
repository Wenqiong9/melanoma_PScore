
# cd ~/Projects/melanoma_pyroptosis/single_cell/RData/
# nohup Rscript 3-3\ -\ used\ -\ computePScore_singleCell_ST_server.R &
# jobs
# [1]+  Running                 nohup Rscript singleCellST_PS_server.R &


# ==============================================================================
# ST
# ==============================================================================
rm(list=ls())
library(BayesSpace)
library(SingleCellExperiment)
library(tidyr)
library(Matrix)
library(GSVA)

load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")
Gene_lists <- list(PScore = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])


# getRDS()  get data
melanoma <- getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")
Data <- counts(melanoma)
data <- logcounts(melanoma) %>% as.matrix()

Gene_lists$PScore[Gene_lists$PScore %in% rownames(data)==FALSE]
rownames(data)[rownames(data) == "PAN3"] <- "NLRP6"

data_E <- GSVA::gsva(data,Gene_lists,method="ssgsea")
data_ES <- t(data_E) %>% data.frame()
data_ES_scale <- scale(data_ES) 
head(data_ES)

save(data_ES,file="2018_thrane_melanoma_ST_PyropScore.RData")

