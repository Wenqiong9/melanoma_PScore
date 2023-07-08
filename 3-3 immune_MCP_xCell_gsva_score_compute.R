

rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/03.immuneMechanism/"))
load("gset_all.add_PScore.Rdata")
# devtools::install_github("IOBR/IOBR")  # https://github.com/IOBR/IOBR
library(IOBR)
source("d:/Work/Codes/CodeLib/3.DataScience/Scoring/gsva_score.R")


# ============================================================================ #
# 1.data prepare ####
# ============================================================================ #

## check NA ##
for(i in 1:7){
  print(sum(is.na(gset_all$gset_data_matrix_symbol[[i]])))
}
## delete NA ##
tmp <- gset_all$gset_data_matrix_symbol[[1]]
tmp <- tmp[!rownames(tmp)%in% c("LINC00882","LINC01716"),]   ## inhouse_data 13 NA
gset_all$gset_data_matrix_symbol[[1]] <- tmp

## require info which EXPRs are array ##
gset_all$GEO_ID
# [1] "inhouse_data"    "TCGA_metastatic" "GSE19234"        "GSE35640"        "GSE54467"        "GSE65904"        "PRJEB23709"   
# RNA-seq, inhouse data, TCGA, PRJEB23709
# ARRAY, GSE35640, GSE19234, GSE54467, GSE65904
array_parameter <- c(rep(FALSE, 2), rep(TRUE, 4), rep(FALSE, 1))


# ============================================================================ #
# 2. "mcpcounter"####
# ============================================================================ #
immune_method <- "mcpcounter"

score_res <- list()
for (i in 1:7){
  print(paste0("Running for dataset ", gset_all$GEO_ID[i], "..."))
  input_exprs <- gset_all$gset_data_matrix_symbol[[i]]
  input_df <- as.data.frame(input_exprs)
 
  immune_score_res <- deconvo_tme(eset = input_df, method = immune_method, arrays = array_parameter[i], perm = 200)
  colnames(immune_score_res) <- gsub("-", "_", gsub("_CIBERSORT", "", colnames(immune_score_res)))
  score_res[[gset_all$GEO_ID[i]]] <- immune_score_res
}

save(score_res, file = "score_res_mcp.RData")

