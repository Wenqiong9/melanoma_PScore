# 
# data prepare
# 
# "inhouse_data" "TCGA_metastatic" "GSE35640" "GSE19234" "GSE54467" "GSE65904" "PRJEB23709"  
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/03.immuneMechanism/"))

library(tibble)
library(ggpubr)
library(magrittr)


# ============================================================================ #
# load data ####
# ============================================================================ #

# gene list
# load(paste0(path, "Data/Gene_list/PRGs_74.RData"))
load((paste0(path, "Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")))
genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])
interest_genes <- genelist$PRGs_prot

# TCGA
mela_exper <- read.table("d:/Work/Data/TCGA_preliminary/SKCM/melanoma_mRNA.txt")  #FPKM
mela_exper <- mela_exper[!(duplicated(mela_exper$symbol)|mela_exper$symbol=="?"),] 
rownames(mela_exper) <- mela_exper$symbol
mela_exper <- mela_exper[-c(1,2)]
mela_exper <- log2(mela_exper+1)

# GSE65904, GSE54467 (metastatic)
load("d:/Projects/Melanoma_Pyroptosis/Melanoma_Pyroptosis_3_GEOvalidation/1.1Melanoma_6Datasets.Rdata")  # 65904_19234_54467_22155

# inhousedata, PRJEB23709, GSE35640 (metastatic)
inhouse_expr <- openxlsx::read.xlsx("d:/Work/Data/inhouse_data(xiangya)/Supplementary tables_HEYi_Ferroptosis.xlsx",sheet=9,startRow = 2,rowNames = TRUE)
ICB_CIBERSORT_PScore <- readr::read_rds(paste0(path, "Results/05.ICB survival and response/ICB_bulk_survival_response/ICB_CIBERSORT_PScore_5_methods.rds.gz"))
load("d:/Work/Data/dara_ICB_secondary/data_exprs_ICB_list_30.RData")

load("d:/Work/Data/dara_ICB_secondary/data_clinical_ICB_list_30.RData")


# ============================================================================ #
# data prepare####
# ============================================================================ #
gset_all <- Melanoma_6Datasets[1:3,]  # GSE65904 GSE19234 GSE54467

# GSE19234 stage III
GSE19234 <- gset_all$gset_data_matrix_symbol[[2]]
GSE19234_pdata <- gset_all$gset_data_pData[[2]]
GSE19234 <- GSE19234[,GSE19234_pdata$geo_accession[GSE19234_pdata$`stage at reccurence:ch1` != "IV"]]
GSE19234_pdata <- GSE19234_pdata[GSE19234_pdata$`stage at reccurence:ch1`!= "IV",]
gset_all$gset_data_matrix_symbol[[2]] <- GSE19234
gset_all$gset_data_pData[[2]] <- GSE19234_pdata

# GSE65904 METASTATIC
GSE65904 <- gset_all$gset_data_matrix_symbol[[1]]
GSE65904_pdata <- gset_all$gset_data_pData[[1]]
GSE65904 <- GSE65904[,GSE65904_pdata$geo_accession[GSE65904_pdata$`tumor stage:ch1` != "Primary"]]
GSE65904_pdata <- GSE65904_pdata[GSE65904_pdata$`tumor stage:ch1` != "Primary",]
gset_all$gset_data_matrix_symbol[[1]] <- GSE65904
gset_all$gset_data_pData[[1]] <- GSE65904_pdata


gset_all <- add_row(gset_all, GEO_ID = c("TCGA_metastatic","inhouse_data","PRJEB23709","GSE35640"))
gset_all <- gset_all[c(5,4,2,7,3,1,6),]
gset_all$GEO_ID
# [1] "inhouse_data"    "TCGA_metastatic"    "GSE19234"    "GSE35640"     "GSE54467"     "GSE65904"     "PRJEB23709"     

# in-house data
load("d:/Projects/Melanoma_Pyroptosis_new/21.inhouseData_GSE120575_bulk30_PScore&response_补充/RData/inhouse_meta.RData")
gset_all$gset_data_matrix_symbol[[1]] <- inhouse_expr
gset_all$gset_data_pData[[1]] <- inhouse_meta

# TCGA
clinic <- read.delim("d:/Work/Data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
melanoma_pdata <- clinic[clinic$type%in%"SKCM",]  # 筛出SKCM数据
gset_all$gset_data_matrix_symbol[[2]] <- mela_exper[,substr(colnames(mela_exper),14,15)=="06"]
melanoma_pdata <- melanoma_pdata[melanoma_pdata$bcr_patient_barcode %in% gsub("\\.","-",substr(colnames(mela_exper),1,12)),]
gset_all$gset_data_pData[[2]] <- melanoma_pdata
 
# GSE35640
gset_all$gset_data_matrix_symbol[[4]] <- data_exprs_ICB_list_30$GSE35640
gset_all$gset_data_pData[[4]] <- data_clinical_ICB_list_30$GSE35640

# PRJEB23709 PRE
PRJEB23709 <-  ICB_CIBERSORT_PScore$gset_data_matrix_symbol[[8]]  
PRJEB23709_pdata <- ICB_CIBERSORT_PScore$gset_data_pData[[8]]  
PRJEB23709 <- PRJEB23709[PRJEB23709_pdata$PrePost == "PRE"]
PRJEB23709_pdata <- PRJEB23709_pdata[PRJEB23709_pdata$SRR %in% colnames(PRJEB23709),]
gset_all$gset_data_matrix_symbol[[7]] <-  PRJEB23709
gset_all$gset_data_pData[[7]] <-  PRJEB23709_pdata



# ============================================================================ #
# check genes ####
# ============================================================================ #

for(i in 1:7){
  ExprMatrix <- gset_all$gset_data_matrix_symbol[[i]]
  rownames(ExprMatrix)[rownames(ExprMatrix)=="PJVK"] <- "DFNB59"
  gset_all$gset_data_matrix_symbol[[i]] <-ExprMatrix
}

ExprMatrix <- gset_all$gset_data_matrix_symbol[[5]]  # GSE54467  
rownames(ExprMatrix)[rownames(ExprMatrix)=="GSDML"] <- "GSDMB"
rownames(ExprMatrix)[rownames(ExprMatrix)=="GSDMDC1"] <- "GSDMD"
rownames(ExprMatrix)[rownames(ExprMatrix)=="CARD12"] <- "NLRC4"
rownames(ExprMatrix)[rownames(ExprMatrix)=="PAN3"] <- "NLRP6"
rownames(ExprMatrix)[rownames(ExprMatrix)=="NALP7"] <- "NLRP7"
rownames(ExprMatrix)[rownames(ExprMatrix)=="CARD15"] <- "NOD2"
gset_all$gset_data_matrix_symbol[[5]] <-ExprMatrix


for(i in 1:7){
  print(interest_genes[interest_genes %in% rownames(gset_all$gset_data_matrix_symbol[[i]])==FALSE])
}
# character(0)
# character(0)
# [1] "CHMP4A"
# [1] "CHMP4A"
# [1] "DFNB59"
# character(0)
# character(0)


## cor - using data after log2
lapply(gset_all$gset_data_matrix_symbol, function(x)head(x)[,1:6])
lapply(gset_all$gset_data_matrix_symbol, function(x)max(x,na.rm=T))

gset_all$gset_data_matrix_symbol[[3]] <- log2(gset_all$gset_data_matrix_symbol[[3]] + 1)
gset_all$gset_data_matrix_symbol[[6]] <- log2(gset_all$gset_data_matrix_symbol[[6]] + 1)
gset_all$gset_data_matrix_symbol[[7]] <- log2(gset_all$gset_data_matrix_symbol[[7]] + 1)

save(gset_all, file = "d:/Work/Projects/BMC_Revised/Results/03.immuneMechanism/gset_all.7_datasets_used.RData")

