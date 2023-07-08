# 
# immune mechanism PScore ####
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/03.immuneMechanism/"))

load("gset_all.7_datasets_used.RData")
load((paste0(path, "Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")))
interest_genes <- uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"]
genelist <- list("PScore" = interest_genes)


# save RData
for(i in 1:7){
  data <- gset_all$gset_data_matrix_symbol[[i]] %>% as.matrix()
  data_ES <- GSVA::gsva(data,genelist,method="ssgsea")
  data_ES <- t(data_ES) %>% data.frame()
  data_ES$PScore_zscore <- scale(data_ES$PScore)
  gset_all[["PScore"]][[i]] <- data_ES
}
save(gset_all,file = "gset_all.add_PScore.Rdata")



# ============================================================================ #
# save xlsx for Table S4
# ============================================================================ #
data_ES_list <- list()
for(i in 1:7){
  data <- gset_all$gset_data_matrix_symbol[[i]] %>% as.matrix()
  data_ES <- GSVA::gsva(data,genelist,method="ssgsea")
  data_ES <- as.data.frame(t(data_ES))
  data_ES$Patient_id <- rownames(data_ES) ###
  data_ES$Dataset <- gset_all$GEO_ID[[i]] ###
  data_ES_list[[gset_all$GEO_ID[[i]]]] <- data_ES
}
openxlsx::write.xlsx(data_ES_list, file = "d:/Work/Projects/BMC_Revised/Results/07.tables/Table S4 PScores of patients.xlsx", overwrite = T)


# load single cell score and spatial score
singleCell_PScore <- readRDS("d:/Work/Projects/BMC_Revised/Results/04.singleCell_ST/singleCell_PScore_GSE115978_GSE75260.rds")
load("d:/Work/Projects/BMC_Revised/Results/04.singleCell_ST/data_ES_spatial.RData")

tmp <- as.data.frame(singleCell_PScore$Pyroptosis_Score[[1]])
tmp$Cell_id <- rownames(tmp)
tmp$Dataset <- "GSE115978"
data_ES_list[["GSE115978"]] <- tmp

tmp <- as.data.frame(singleCell_PScore$Pyroptosis_Score[[2]])
tmp$Cell_id <- rownames(tmp)
tmp$Dataset <- "GSE72056"
data_ES_list[["GSE72056"]] <- tmp

data_ES$cell_id <- rownames(data_ES)
data_ES$Dataset <- "2018_thrane_melanoma"
data_ES_list[["2018_thrane_melanoma"]] <- data_ES

save(data_ES_list, file = "d:/Work/Projects/BMC_Revised/Results/07.tables/Table S4 PScores of patients.RData")
openxlsx::write.xlsx(data_ES_list, file = "d:/Work/Projects/BMC_Revised/Results/07.tables/Table S4 PScores of patients.xlsx", overwrite = T)
