# cd ~/Projects/melanoma_pyroptosis/single_cell/RData/
# nohup Rscript 3-3\ -\ used\ -\ computePScore_singleCell_ST_server.R &
# jobs
# [1]+  Running                 nohup Rscript singleCellST_PS_server.R &


# ==============================================================================
# compute score ####
# ==============================================================================
rm(list=ls())
library(Seurat)
library(harmony)
library(Rcpp)
library(tidyr)

# SKCM_TISHdata <- readRDS('/home/pub/Data/public_scRNA/TISH_DB_Melanoma/TISH_DB_ScRNAseqData.rds.gz')
SKCM_TISHdata <- readRDS('d:/Work/Data/public_scRNA/TISH_DB_Melanoma/TISH_DB_ScRNAseqData.rds.gz')
SKCM_TISHdata <- SKCM_TISHdata[c(1,6),]

load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")
interest_genes <- uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"]
Gene_lists <- list(PScore=interest_genes)


MyGenes <- Gene_lists$PScore
for(i in 1:2){
  print(paste0("Running for Task",i))
  SeuratData <- SKCM_TISHdata$SeuratObject[[i]]
  rownames(SeuratData@assays$RNA@data)[rownames(SeuratData@assays$RNA@data)=="PJVK"] <- "DFNB59"
  print(table(MyGenes %in% rownames(SeuratData@assays$RNA@data)))
  print(MyGenes[MyGenes %in% rownames(SeuratData@assays$RNA@data)==FALSE])
}


# GSE115978   29 metas + 4 prima;  pre post treatment
# memory.limit(110000)

Pyroptosis_Score_melanoma <- SKCM_TISHdata %>% dplyr::mutate(Pyroptosis_Score=purrr::map(.x=SeuratObject,function(.x){
  SeuratData <- .x
  matrix <- SeuratData@assays$RNA@data %>% as.matrix()
  data_ES <- GSVA::gsva(matrix,Gene_lists,method="ssgsea")
  data_ES <- t(data_ES) %>% data.frame()
  return(data_ES)
})) %>% dplyr::select(-SeuratObject)

# Save Data
readr::write_rds(Pyroptosis_Score_melanoma,file="d:/Work/Projects/BMC_Revised/Results/04.singleCell_ST/singleCell_PScore_GSE115978_GSE75260.rds")
