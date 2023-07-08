# 
# aim: multi methods -- compare ROC of survival
# 
# ============================================================================ #

rm(list=ls())
options(stringsAsFactors = F)
library(GSVA)
library(dplyr)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
melanoma <- read.table("d:/Projects/Melanoma_Pyroptosis/Genelist/melanoma_mRNA.txt")  #FPKM
load("coxForest/uniCoxSig_list.RData")


# ============================================================================ #
# computing PScore ####
# ============================================================================ #
# prepare matrix
melanoma[1:4,1:4]
mela_exper <- melanoma[melanoma$symbol!="?",] 
mela_exper <- mela_exper[!duplicated(mela_exper$symbol),]
rownames(mela_exper) <- mela_exper$symbol
mela_exper <- mela_exper[-c(1,2)]
head(mela_exper)[,1:6]; dim(mela_exper)
mela_exper <- log2(mela_exper+1)
mela_exper <- mela_exper[,!substr(colnames(mela_exper),14,15) > 10] 

# PScore -- metastatic -- ssgsea
melanoma_expr_meta <- mela_exper[,which(substr(colnames(mela_exper),14,15) == "06")]
genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])
genelist$PRGs_prot[genelist$PRGs_prot %in% rownames(melanoma_expr_meta) ==FALSE]

data_ES_met <- GSVA::gsva(as.matrix(melanoma_expr_meta),genelist,method="ssgsea")
data_ES_met_t <- as.data.frame(t(data_ES_met))
head(data_ES_met_t)
data_ES_met_t_zscore <- scale(data_ES_met_t) %>% data.frame()

save(data_ES_List,data_ES_zscore_List,file = "PScoreModel_ROC/PyropScore_SKCM.Rdata")