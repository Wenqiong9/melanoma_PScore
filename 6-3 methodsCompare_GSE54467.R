# GEO median survival analysis & multiCOX
# Methods schematic. Tissue samples from patients with AJCC stage III melanoma 
# were obtained from a prospective collection of fresh-frozen tumors. 
# https://onlinelibrary.wiley.com/doi/10.1002/ijc.29047
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
library(survival)
library(dplyr)
library(GSVA)
# path <- "d:/Work/Projects/BMC_Revised/"
# setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
load("d:/Work/Data/GEO/GEO_melanoma_res/Melanoma_Datasets_16.RData")
load("d:/Work/Projects/BMC_Revised/Data/Gene_list/PRGs_74.RData")
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")


# ============================================================================ #
# prepare ####
# ============================================================================ #
pData <- Melanoma_Datasets_16$gset_data_pData[[11]]
pData <- pData[,c(2,(ncol(pData)-6):ncol(pData))]

pData$OS.time <- pData$`survival from primary melanoma (months):ch1` %>% as.numeric()
pData$OS_year <- as.numeric(pData$OS.time)/12

pData$OS.status <- pData$`patient last status:ch1`
pData$OS.status[which(pData$OS.status %in% "Dead Melanoma")]<-1
pData$OS.status[which(pData$OS.status %in% c("Alive NSR","Alive with Melanoma"))]<-0
pData$OS.status[which(pData$OS.status %in% c("Alive Status Unknown","Dead Cause Unknown","Dead Not Melanoma"))]<-NA
pData$OS.status <- as.numeric(pData$OS.status)

ExprMatrix <- Melanoma_Datasets_16$gset_data_matrix_symbol[[11]]


# ============================================================================ #
# GSE54467: PScore ####
# ============================================================================ #
genelist <- list(PRGs_74 = PRGs_74$Symbol.ID,
                 PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"],
                 Wang_Ann.Transl.Med_2022 = c("AIM2", "GSDMC", "GSDMD", "IL18", "NLRP6", "PRKACA"))

genelist$PRGs_prot[which(genelist$PRGs_prot %in% rownames(ExprMatrix)==FALSE)]
rownames(ExprMatrix)[rownames(ExprMatrix)=="CARD15"] <- "NOD2"
rownames(ExprMatrix)[rownames(ExprMatrix)=="GSDML"] <- "GSDMB"
rownames(ExprMatrix)[rownames(ExprMatrix)=="GSDMDC1"] <- "GSDMD"
rownames(ExprMatrix)[rownames(ExprMatrix)=="CARD12"] <- "NLRC4"
rownames(ExprMatrix)[rownames(ExprMatrix)=="PAN3"] <- "NLRP6"
rownames(ExprMatrix)[rownames(ExprMatrix)=="NALP7"] <- "NLRP7"
rownames(ExprMatrix)[rownames(ExprMatrix)=="MLZE"] <- "GSDMC"
# rownames(ExprMatrix)[rownames(ExprMatrix)=="PJVK"]  # no GSDMF PJVK DFNB59

genelist$PRGs_prot[which(genelist$PRGs_prot %in% rownames(ExprMatrix)==FALSE)]
genelist$Wang_Ann.Transl.Med_2022[which(genelist$Wang_Ann.Transl.Med_2022 %in% rownames(ExprMatrix)==FALSE)]
# "PRKACA"

data_ES <- GSVA::gsva(as.matrix(ExprMatrix),genelist,method="ssgsea")
data_ES <- t(data_ES) %>% data.frame()
data_ES_zscore <- scale(data_ES) %>% data.frame()  


# coefficient
ExprMatrix <- as.data.frame(t(ExprMatrix))

# Zhu_Stem.Cells.Int_2023: Score=(-0.166*CCL8)-0.156*FCGR2A+0.047*GBP2-0.327*GRIPAP1-0.207*HAPLN3+0.145*HPDL-0.022*IFITM1-1.01*TRIM34.
# c("HPDL", "4-HPPD-L", "GLOXD1", "NEDSWMA", "SPG83") %in% colnames(ExprMatrix)
# c("TRIM34", "IFP1", "RNF21") %in% colnames(ExprMatrix)
# ExprMatrix$Zhu_Stem.Cells.Int_2023 <- (-0.166 * ExprMatrix$CCL8) - 0.156 * ExprMatrix$FCGR2A +
#     0.047 * ExprMatrix$GBP2 - 0.327 * ExprMatrix$GRIPAP1 -
#     0.207 * ExprMatrix$HAPLN3 + 0.145 * ExprMatrix$GLOXD1 -  # HPDL: 4-HPPD-L, GLOXD1, NEDSWMA, SPG83
#     0.022 * ExprMatrix$IFITM1 - 1.01 * ExprMatrix$TRIM34   # TRIM34: IFP1, RNF21
ExprMatrix$Zhu_Stem.Cells.Int_2023 <- (-0.166 * ExprMatrix$CCL8) - 0.156 * ExprMatrix$FCGR2A +
    0.047 * ExprMatrix$GBP2 - 0.327 * ExprMatrix$GRIPAP1 -
    0.207 * ExprMatrix$HAPLN3 + 0.145 * ExprMatrix$GLOXD1 - 
    0.022 * ExprMatrix$IFITM1   # no TRIM34

# Li_Int.J.Gen.Med_2022: Score=(-0.119*AIM2)+(-0.487*NLRP6)+(-0.374*IL18)+(0.230*GSDMA)+(0.383*GSDMC)
# c("GSDMA", "FKSG9", "GSDM", "GSDM1") %in% colnames(ExprMatrix) 
ExprMatrix$Li_Int.J.Gen.Med_2022 <- (-0.119 * ExprMatrix$AIM2) + (-0.487 * ExprMatrix$NLRP6)+
    (-0.374 * ExprMatrix$IL18) + (0.230 * ExprMatrix$GSDM1)+
    (0.383 * ExprMatrix$GSDMC)

# Xu_Front.Med_2022: Score=0.003*EMP3-0.065*TLR1-0.012*IFNGR2-0.288*IL15-0.057*CCL8-0.633*NLRP6-0.329*CCL25-0.024*RTP4
ExprMatrix$Xu_Front.Med_2022 <- (0.003 * ExprMatrix$EMP3) - 0.065 * ExprMatrix$TLR1 - 
    0.012 * ExprMatrix$IFNGR2 - 0.288 * ExprMatrix$IL15 - 
    0.057 * ExprMatrix$CCL8 - 0.633 * ExprMatrix$NLRP6 - 
    0.329 * ExprMatrix$CCL25 - 0.024 * ExprMatrix$RTP4

# Niu_Math.Biosci.Eng_2022: Score=0.038*GSDMA+0.31*GSDMC-0.028*AIM2-0.437*NOD2
ExprMatrix$Niu_Math.Biosci.Eng_2022 <- 0.038 * ExprMatrix$GSDM1 + 0.31 * ExprMatrix$GSDMC - 
    0.028 * ExprMatrix$AIM2 - 0.437 * ExprMatrix$NOD2

# Wu_PeerJ_2021: Score=0.2758*GSDMC-0.0699*GZMA-0.0526*AIM2-0.1766*PDL1
ExprMatrix$Wu_PeerJ_2021 <- 0.2758 * ExprMatrix$GSDMC - 0.0699 * ExprMatrix$GZMA - 
    0.0526 * ExprMatrix$AIM2 - 0.1766 * ExprMatrix$CD274
# Ju_Front.Oncol_2021: Score=-0.006861*GSDMD+0.0003969*GSDME-0.001943*CASP4+0.0079361*GSDMC-0.022123*NLRC4-0.009636*APIP-0.003569*AIM2-0.00106*CASP3-0.000169*IL18
ExprMatrix$Ju_Front.Oncol_2021 <- -0.006861 * ExprMatrix$GSDMD + 0.0003969 * ExprMatrix$DFNA5 - 
    0.001943 * ExprMatrix$CASP4 + 0.0079361 * ExprMatrix$GSDMC - 
    0.022123 * ExprMatrix$NLRC4 - 0.009636 * ExprMatrix$APIP - 
    0.003569 * ExprMatrix$AIM2 - 0.00106 * ExprMatrix$CASP3 - 
    0.000169 * ExprMatrix$IL18

# Shi_Medicine_2022: Score=-0.0452*BST2-0.1636*GBP5-0.0531*AIM2
ExprMatrix$Shi_Medicine_2022 <-  -0.0452 * ExprMatrix$BST2 - 0.1636 * ExprMatrix$GBP5 - 0.0531 * ExprMatrix$AIM2

# Wu_Cancer.Med_2022: Score=0.139*CASP5+0.240*NLRP6+1.388*NLRP7+0.112*PYCARD
ExprMatrix$Wu_Cancer.Med_2022 <- 0.139 * ExprMatrix$CASP5 + 0.240 * ExprMatrix$NLRP6 + 
    1.388 * ExprMatrix$NLRP7 + 0.112 * ExprMatrix$PYCARD

# Wang_J.Investig.Dermatol_2022: Score = mean(IRF9, STAT2)
# c("IRF-9", "ISGF3", "ISGF3G", "p48") %in% colnames(ExprMatrix)
ExprMatrix$Wang_J.Investig.Dermatol_2022 <- apply(ExprMatrix[,c("ISGF3G","STAT2")],1,mean)


head(data_ES)
head(ExprMatrix[,20503:20511])
data_ES$ID <- rownames(data_ES)
data_ES$ID == rownames(ExprMatrix)

multiScore <- cbind(data_ES,ExprMatrix[,20503:20511])
multiScore$ID == pData$geo_accession
GSE54467_m <- merge(multiScore,pData,by.x="ID",by.y="geo_accession")


# ============================================================================ #
# 2.cox forest ####
# ============================================================================ #
library(survival)
library(ggplot2)

uniCox <- data.frame()
for(j in colnames(GSE54467_m)[2:13]){
    cox <- coxph(Surv(OS.time,OS.status) ~ GSE54467_m[,j], data = GSE54467_m)
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
uniCox$Risk_group <- ifelse(uniCox$HR >= 1 & uniCox$P_value < 0.05, "Risk",
                            ifelse(uniCox$HR < 1 & uniCox$P_value < 0.05, "Protect","Not_sig"))

ggplot(uniCox, aes(HR, gene)) + 
    geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
    geom_point(shape=18, aes(col=Risk_group, size = P_value))  +
    geom_vline(aes(xintercept = 1),color="darkgrey",linetype = "dashed",size=0.8) + 
    scale_x_log10() +
    theme_bw() 

uniCox_GSE54467 <- uniCox
uniCox_GSE54467$dataset <- "GSE54467"
save(uniCox_GSE54467, file = "d:/Work/Projects/BMC_Revised/Results/08.methods_compare/uniCox_GSE54467.RData")
