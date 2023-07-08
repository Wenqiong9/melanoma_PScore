# aim: ICB 


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/05.ICB survival and response/"))
library(dplyr)
ICB_CIBERSORT <- readr::read_rds("d:/Work/Data/public_immune/Datasets/ICB_CIBERSORT_used13.rds.gz")
load(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData"))
genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])

# for(i in c(1:13)){
#   ExprMatrix <- ICB_CIBERSORT$gset_data_matrix_symbol[[i]]
#   print(genelist$PRGs_prot[genelist$PRGs_prot %in% rownames(ExprMatrix) == FALSE])
# }
# ExprMatrix <- ICB_CIBERSORT$gset_data_matrix_symbol[[1]]
# c("C14orf123", "CHMP4", "CHMP4B", "HSPC134", "SHAX2", "SNF7", "SNF7-1", "VPS32-1", "VPS32A")  %in% rownames(ExprMatrix)


# ============================================================================ #
# PRGs_prot ####
# ============================================================================ #
ICB_CIBERSORT <- ICB_CIBERSORT %>% dplyr::mutate(Pyroptosis_score=purrr::map(.x=gset_data_matrix_symbol,function(.x){
  Data <- .x %>% as.matrix()
  data_ES <- GSVA::gsva(Data,genelist,method="ssgsea")
  data_ES <- data_ES %>% t() %>% data.frame()
  data_ES$id <- rownames(data_ES)
  return(data_ES)
}))

ICB_CIBERSORT <- ICB_CIBERSORT %>% dplyr::mutate(Pyroptosis_score_pdata=purrr::map2(.x=Pyroptosis_score,.y=gset_data_pData,function(.x,.y){
  Data <- .x %>% data.frame()
  Data2 <- .y %>% data.frame()
  Data <- merge(Data, Data2, by.x="id", by.y="sample_use")
  return(Data)
}))

# readr::write_rds(ICB_CIBERSORT,file = "ICB_bulk_survival_response/ICB_CIBERSORT_PScore_Sig34.rds.gz",compress = "gz")
readr::write_rds(ICB_CIBERSORT,file = "ICB_bulk_survival_response/ICB_CIBERSORT_PScore.rds.gz",compress = "gz")
