# 
# aim: multi methods -- forest cox
# 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
library(GSVA)
library(dplyr)
melanoma <- read.table("d:/Projects/Melanoma_Pyroptosis/Genelist/melanoma_mRNA.txt")  #FPKM
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")


# ============================================================================ #
# 1. computing PScore ####
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
genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"],
                 Wang_Ann.Transl.Med_2022 = c("AIM2", "GSDMC", "GSDMD", "IL18", "NLRP6", "PRKACA"))
genelist$PRGs_prot[genelist$PRGs_prot %in% rownames(mela_exper) == FALSE]

data_ES_met <- GSVA::gsva(as.matrix(melanoma_expr_meta),genelist,method="ssgsea")
data_ES_met_t <- t(data_ES_met) %>% data.frame()
head(data_ES_met_t)
data_ES_met_t_zscore <- scale(data_ES_met_t) %>% data.frame()


# coefficient
melanoma_expr_meta <- mela_exper[,which(substr(colnames(mela_exper),14,15) == "06")]
melanoma_expr_meta <- t(melanoma_expr_meta) %>% data.frame()
# Zhu_Stem.Cells.Int_2023: Score=(-0.166*CCL8)-0.156*FCGR2A+0.047*GBP2-0.327*GRIPAP1-0.207*HAPLN3+0.145*HPDL-0.022*IFITM1-1.01*TRIM34.
melanoma_expr_meta$Zhu_Stem.Cells.Int_2023 <- (-0.166 * melanoma_expr_meta$CCL8) - 0.156 * melanoma_expr_meta$FCGR2A + 
  0.047 * melanoma_expr_meta$GBP2 - 0.327 * melanoma_expr_meta$GRIPAP1 - 
  0.207 * melanoma_expr_meta$HAPLN3 + 0.145 * melanoma_expr_meta$HPDL - 
  0.022 * melanoma_expr_meta$IFITM1 - 1.01 * melanoma_expr_meta$TRIM34
# Li_Int.J.Gen.Med_2022: Score=(-0.119*AIM2)+(-0.487*NLRP6)+(-0.374*IL18)+(0.230*GSDMA)+(0.383*GSDMC)
melanoma_expr_meta$Li_Int.J.Gen.Med_2022 <- (-0.119 * melanoma_expr_meta$AIM2) + (-0.487 * melanoma_expr_meta$NLRP6)+
  (-0.374 * melanoma_expr_meta$IL18) + (0.230 * melanoma_expr_meta$GSDMA)+
  (0.383 * melanoma_expr_meta$GSDMC)
# Xu_Front.Med_2022: Score=0.003*EMP3-0.065*TLR1-0.012*IFNGR2-0.288*IL15-0.057*CCL8-0.633*NLRP6-0.329*CCL25-0.024*RTP4
melanoma_expr_meta$Xu_Front.Med_2022 <- (0.003 * melanoma_expr_meta$EMP3) - 0.065 * melanoma_expr_meta$TLR1 - 
  0.012 * melanoma_expr_meta$IFNGR2 - 0.288 * melanoma_expr_meta$IL15 - 
  0.057 * melanoma_expr_meta$CCL8 - 0.633 * melanoma_expr_meta$NLRP6 - 
  0.329 * melanoma_expr_meta$CCL25 - 0.024 * melanoma_expr_meta$RTP4
# Niu_Math.Biosci.Eng_2022: Score=0.038*GSDMA+0.31*GSDMC-0.028*AIM2-0.437*NOD2
melanoma_expr_meta$Niu_Math.Biosci.Eng_2022 <- 0.038 * melanoma_expr_meta$GSDMA + 0.31 * melanoma_expr_meta$GSDMC - 
  0.028 * melanoma_expr_meta$AIM2 - 0.437 * melanoma_expr_meta$NOD2
# Wu_PeerJ_2021: Score=0.2758*GSDMC-0.0699*GZMA-0.0526*AIM2-0.1766*PDL1
melanoma_expr_meta$Wu_PeerJ_2021 <- 0.2758 * melanoma_expr_meta$GSDMC - 0.0699 * melanoma_expr_meta$GZMA - 
  0.0526 * melanoma_expr_meta$AIM2 - 0.1766 * melanoma_expr_meta$CD274
# Ju_Front.Oncol_2021: Score=-0.006861*GSDMD+0.0003969*GSDME-0.001943*CASP4+0.0079361*GSDMC-0.022123*NLRC4-0.009636*APIP-0.003569*AIM2-0.00106*CASP3-0.000169*IL18
melanoma_expr_meta$Ju_Front.Oncol_2021 <- -0.006861 * melanoma_expr_meta$GSDMD + 0.0003969 * melanoma_expr_meta$DFNA5 - 
  0.001943 * melanoma_expr_meta$CASP4 + 0.0079361 * melanoma_expr_meta$GSDMC - 
  0.022123 * melanoma_expr_meta$NLRC4 - 0.009636 * melanoma_expr_meta$APIP - 
  0.003569 * melanoma_expr_meta$AIM2 - 0.00106 * melanoma_expr_meta$CASP3 - 
  0.000169 * melanoma_expr_meta$IL18
# Shi_Medicine_2022: Score=-0.0452*BST2-0.1636*GBP5-0.0531*AIM2
melanoma_expr_meta$Shi_Medicine_2022 <-  -0.0452 * melanoma_expr_meta$BST2 - 0.1636 * melanoma_expr_meta$GBP5 - 0.0531 * melanoma_expr_meta$AIM2
# Wu_Cancer.Med_2022: Score=0.139*CASP5+0.240*NLRP6+1.388*NLRP7+0.112*PYCARD
melanoma_expr_meta$Wu_Cancer.Med_2022 <- 0.139 * melanoma_expr_meta$CASP5 + 0.240 * melanoma_expr_meta$NLRP6 + 
  1.388 * melanoma_expr_meta$NLRP7 + 0.112 * melanoma_expr_meta$PYCARD
# Wang_J.Investig.Dermatol_2022: Score = mean(IRF9, STAT2)
melanoma_expr_meta$Wang_J.Investig.Dermatol_2022 <- apply(melanoma_expr_meta[,c("IRF9","STAT2")],1,mean)


head(data_ES_met_t)
head(melanoma_expr_meta[,20502:20510])

multiScore <- cbind(data_ES_met_t,melanoma_expr_meta[,20502:20510])


# ============================================================================ #
# 2.cox -- metastatic####
# ============================================================================ #
library(survival)
library(ggplot2)
load("D:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/Expr_Survival_list.PRG_74.RData")

melanoma_clinic <- Expr_Survival_list$metastatic
melanoma_clinic$sample_id <- rownames(melanoma_clinic)
multiScore$sample_id <- gsub("\\.","-",substr(rownames(multiScore),1,12))
multiScore_clinic_merge <- merge(multiScore, melanoma_clinic[,c(1,2,77)], by = "sample_id")
multiScore_clinic_merge$OS_year <- multiScore_clinic_merge$OS.time/365

uniCox <- data.frame()
# j = colnames(multiScore_clinic_merge)[3]
# j = colnames(multiScore_clinic_merge)[4]
# j = colnames(multiScore_clinic_merge)[11]
for(j in colnames(multiScore_clinic_merge)[2:12]){
  cox <- coxph(Surv(OS.time,OS) ~ multiScore_clinic_merge[,j], data = multiScore_clinic_merge)
  coxSummary <- summary(cox)
  C_index = signif(coxSummary$concordance[["C"]],3)
  uniCox <- rbind(uniCox, data.frame(gene = j,
                                     HR = coxSummary$conf.int[,"exp(coef)"],
                                     HR.95L = coxSummary$conf.int[,"lower .95"],
                                     HR.95H = coxSummary$conf.int[,"upper .95"],
                                     P_value = coxSummary$coefficients[,"Pr(>|z|)"],
                                     Concordance = C_index,
                                     Concordance_se = paste0(C_index, "(se = ", signif(coxSummary$concordance[["se(C)"]],3),")")))
}
uniCox$Risk_group <- ifelse(uniCox$HR > 1 & uniCox$P_value < 0.05, "Risk",
                            ifelse(uniCox$HR < 1 & uniCox$P_value < 0.05, "Protect","Not_sig"))

uniCox_TCGA <- uniCox
uniCox_TCGA$dataset <- "TCGA_metastatic"
save(uniCox_TCGA, file = "d:/Work/Projects/BMC_Revised/Results/08.methods_compare/uniCox_TCGA.RData")

ggplot(uniCox, aes(HR, gene)) + 
  geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
  geom_point(shape=18, aes(col=Risk_group, size = P_value))  +
  geom_vline(aes(xintercept = 1),color="darkgrey",linetype = "dashed",size=0.8) + 
  scale_x_log10() +
  theme_bw() 
