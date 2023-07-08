# Aim ：data prepare -- single gene COX -- HR forest
# expression matrix: fpkm
# ============================================================================ #


rm(list=ls())
options("scipen"=100,"digits"=7)
library(dplyr)
library(survival)
library(ggplot2)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))

mela_exper <- read.table("d:/Projects/Melanoma_Pyroptosis/Genelist/melanoma_mRNA.txt")  #FPKM
clinic <- read.delim("d:/Projects/Melanoma_Pyroptosis/Genelist/TCGA_ClinicalData_20180420.txt")
load(paste0(path, "Data/Gene_list/PRGs_74.RData"))


# ============================================================================ #
# 1. data_prepare: merge gene Expr and melanoma clinic Info ####
# ============================================================================ #

## 1.1 Extract the expression matrix of interested genes ####

interest_genes <- PRGs_74
melanoma_interest <- mela_exper[mela_exper$symbol%in% c(interest_genes$Symbol.ID),]
rownames(melanoma_interest) <- melanoma_interest$symbol
melanoma_interest <- melanoma_interest[-c(1,2)]
head(melanoma_interest)[,1:6]

# log2 transform
melanoma_interest <- log2(melanoma_interest+1)

# Extract the expression matrix of patients with primary and normal tumors 
# Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29.
melanoma_interest <- melanoma_interest[,!substr(colnames(melanoma_interest),14,15) > 10] 
melanoma_expr_list <- list(tumor = melanoma_interest,
                           primary = melanoma_interest[,which(substr(colnames(melanoma_interest),14,15) == "01")],
                           metastatic = melanoma_interest[,which(substr(colnames(melanoma_interest),14,15) == "06")])

save(melanoma_expr_list, file = "coxForest/melanoma_expr_list.PRGs_74.RData")


## 1.2 merge clinical infomation - OS, OS.status ####

SKCM_clinic <- clinic[clinic$type %in% "SKCM",c("bcr_patient_barcode","OS","OS.time")]
colnames(SKCM_clinic) <- c("id","OS","OS.time")

Expr_Survival_list <- list()
for(i in c("primary", "metastatic")){
  Expr_tmp <- melanoma_expr_list[[i]]
  
  Expr_tmp <- t(Expr_tmp) %>% data.frame()
  Expr_tmp$id <- gsub("\\.", "-", substr(rownames(Expr_tmp), 1, 12))
  Expr_tmp <- Expr_tmp[!duplicated(Expr_tmp$id),]    
  Survtime_Expr <- merge(SKCM_clinic, Expr_tmp)     
  Survtime_Expr[,2:3] <- apply(Survtime_Expr[,2:3],2,as.numeric)
  Survtime_Expr <- Survtime_Expr[is.na(Survtime_Expr$OS)==FALSE,]  # patient whose OS are NA also without OS.time
  rownames(Survtime_Expr) <- Survtime_Expr$id
  Survtime_Expr <- Survtime_Expr[, !colnames(Survtime_Expr) %in% "id"]

  Expr_Survival_list[[i]] <- Survtime_Expr
}

save(Expr_Survival_list,file = "coxForest/Expr_Survival_list.PRG_74.RData")


# ============================================================================ #
# 2. screening for prognostic genes by single factor COX regression ####
# ============================================================================ #

uniCox_list <- list()
uniCoxSig_list <- list()
uniCoxSig_Exp_list <- list() 
for(i in names(Expr_Survival_list)){
  print(paste0("Running for Task ",i))
  Survtime_Expr <- Expr_Survival_list[[i]]
  colnames(Survtime_Expr)[c(1,2)] <- c("surv.status","surv.time")

  # cycle - single gene cox
  pFilter <- 0.05 
  uniCox <- data.frame()
  uniCoxSig <- data.frame()
  for(j in colnames(Survtime_Expr[,3:ncol(Survtime_Expr)])){
    cox <- coxph(Surv(surv.time,surv.status) ~ Survtime_Expr[,j], data = Survtime_Expr)
    coxSummary <- summary(cox)
    uniCox <- rbind(uniCox, data.frame(gene = j,
                                       HR = coxSummary$conf.int[,"exp(coef)"],
                                       HR.95L = coxSummary$conf.int[,"lower .95"],
                                       HR.95H = coxSummary$conf.int[,"upper .95"],
                                       P_value = coxSummary$coefficients[,"Pr(>|z|)"]))
    coxP <- coxSummary$coefficients[,"Pr(>|z|)"]
    if(coxP<pFilter){
      uniCoxSig <- rbind(uniCoxSig, data.frame(gene = j, 
                                               HR = coxSummary$conf.int[,"exp(coef)"],
                                               HR.95L = coxSummary$conf.int[,"lower .95"],
                                               HR.95H = coxSummary$conf.int[,"upper .95"],
                                               P_value = coxSummary$coefficients[,"Pr(>|z|)"]))
      }
    }
  
  # add prognostic label
  uniCox$Risk_group <- ifelse(uniCox$HR >= 1 & uniCox$P_value < 0.05, "Risk",
                              ifelse(uniCox$HR < 1 & uniCox$P_value < 0.05, "Protect","Not_sig"))
  uniCoxSig$Risk_group <- ifelse(uniCoxSig$HR >= 1, "Risk", "Protect")
  print(uniCoxSig$gene)
    
  # save matrix of prognostic genes
  uniCoxSig_Exp <- Survtime_Expr[,c("surv.status", "surv.time", uniCoxSig$gene)]
  uniCoxSig_Exp <- cbind(id = row.names(uniCoxSig_Exp), uniCoxSig_Exp)
  
  uniCox_list[[i]] <- uniCox
  uniCoxSig_list[[i]] <- uniCoxSig
  uniCoxSig_Exp_list[[i]] <- uniCoxSig_Exp
}

save(uniCox_list,file = "coxForest/uniCox_list.PRG_74.RData")
save(uniCoxSig_list,file = "coxForest/uniCoxSig_list.RData")
save(uniCoxSig_Exp_list,file = "coxForest/uniCoxSig_Exp_list.RData")

# [1] "Running for Task primary"
# [1] "IRF1" "IRF2" "ZBP1"
# [1] "Running for Task metastatic"
# [1] "AIM2"     "APIP"     "BAK1"     "BAX"      "CASP1"    "CASP3"    "CASP4"    "CASP5"    "CASP8"    "CFLAR"    "CHMP2B"  
# [12] "CHMP4A"   "CHMP5"    "DFNB59"   "DPP9"     "GSDMB"    "GSDMD"    "GZMA"     "GZMB"     "IL18"     "IL1B"     "IRF1"    
# [23] "IRF2"     "MEFV"     "NAIP"     "NLRC4"    "NLRP1"    "NLRP3"    "NLRP6"    "NLRP7"    "NOD2"     "TLR4"     "TNFRSF1B"
# [34] "ZBP1"  


# ============================================================================ #
# 3. forest plot ####
# ============================================================================ #

names(uniCox_list)<- c("SKCM (Primary)","SKCM (Metastastic)")
openxlsx::write.xlsx(uniCox_list, "coxForest/univariateCOX_OS.PRG_74_primary_metastatic.xlsx")

mydata <- uniCox_list$`SKCM (Primary)`
ggplot(mydata, aes(HR, gene)) + 
  geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
  geom_point(shape=18, aes(col=Risk_group, size = p_value))  +
  geom_vline(aes(xintercept = 1),color="darkgrey",linetype = "dashed",size=0.8) + 
  theme_bw() 


