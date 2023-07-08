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
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
load("d:/Work/Data/GEO/GEO_melanoma_res/Melanoma_Datasets_16.RData")
load("coxForest/uniCoxSig_list.RData")
genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])


# ============================================================================ #
# prepare ####
# ============================================================================ #
pData <- Melanoma_Datasets_16$gset_data_pData[[11]]
pData <- pData[,c(2,(ncol(pData)-6):ncol(pData))]

pData$OS.time <- pData$`survival from primary melanoma (months):ch1` %>% as.numeric()
pData$OS.time.year <- as.numeric(pData$OS.time)/12

pData$OS.status <- pData$`patient last status:ch1`
pData$OS.status[which(pData$OS.status %in% "Dead Melanoma")]<-1
pData$OS.status[which(pData$OS.status %in% c("Alive NSR","Alive with Melanoma"))]<-0
pData$OS.status[which(pData$OS.status %in% c("Alive Status Unknown","Dead Cause Unknown","Dead Not Melanoma"))]<-NA
pData$OS.status <- as.numeric(pData$OS.status)

ExprMatrix <- Melanoma_Datasets_16$gset_data_matrix_symbol[[11]]


# ============================================================================ #
# GSE54467: PScore ####
# ============================================================================ #
genelist$PRGs_prot[which(genelist$PRGs_prot %in% rownames(ExprMatrix)==FALSE)]
# [1] "DFNB59" "GSDMB"  "GSDMD"  "NLRC4"  "NLRP6"  "NLRP7"  "NOD2"  
rownames(ExprMatrix)[rownames(ExprMatrix)=="CARD15"] <- "NOD2"
rownames(ExprMatrix)[rownames(ExprMatrix)=="GSDML"] <- "GSDMB"
rownames(ExprMatrix)[rownames(ExprMatrix)=="GSDMDC1"] <- "GSDMD"
rownames(ExprMatrix)[rownames(ExprMatrix)=="CARD12"] <- "NLRC4"
rownames(ExprMatrix)[rownames(ExprMatrix)=="PAN3"] <- "NLRP6"
rownames(ExprMatrix)[rownames(ExprMatrix)=="NALP7"] <- "NLRP7"
# rownames(ExprMatrix)[rownames(ExprMatrix)=="PJVK"]  # no GSDMF PJVK DFNB59

data_ES <- GSVA::gsva(as.matrix(ExprMatrix),genelist,method="ssgsea")
data_ES <- t(data_ES) %>% data.frame()
data_ES_zscore <- scale(data_ES) %>% data.frame()  

data_ES$ID <- rownames(data_ES)
GSE54467_m <- merge(data_ES,pData,by.x="ID",by.y="geo_accession")


# ============================================================================ #
# GSE54467: survival analysis ####
# ============================================================================ #
library(survminer)
library(survival)

pdataScoreMerge <- GSE54467_m

res.cut <- surv_cutpoint(pdataScoreMerge, time = "OS.time", event = "OS.status",variables = "PRGs_prot",minprop = 0.4)
res.cat <- surv_categorize(res.cut)
pdataScoreMerge$Group <- ifelse(res.cat$PRGs_prot == "high","High_PS","Low_PS")
pdataScoreMerge$Group <- factor(pdataScoreMerge$Group,levels = c("Low_PS","High_PS"))

diff=survdiff(Surv(OS.time,OS.status) ~ Group,data = pdataScoreMerge)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue=format(pValue, scientific = TRUE)
pValue

cox.res=coxph(Surv(OS.time,OS.status)~ Group,data = pdataScoreMerge)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)
CI_low <- signif(CI[1],2)
CI_high <- signif(CI[2],2)

fit <- survfit(Surv(OS.time.year,OS.status) ~ Group, data = pdataScoreMerge)
    
ggsurvplot(fit, data=pdataScoreMerge,
           conf.int=F,
           pval=paste0("p=",pValue,"\nHR=",paste0(HR,"(",CI_low,"-",CI_high,")")),
           risk.table=TRUE,
           legend.labs=factor(c("Low_PS","High_PS")),
           legend.title="")


# ============================================================================ #
# GSE54467: Multivariate cox regression ####
# ============================================================================ #
pdataScoreMerge <- pdataScoreMerge[which(pdataScoreMerge$OS.time!="" & pdataScoreMerge$OS.status!="" ),] 
pdataScoreMerge$Gender <- pdataScoreMerge$`patient sex:ch1`
pdataScoreMerge$Gender <- factor(pdataScoreMerge$Gender,levels = c("Male","Female"))
pdataScoreMerge$Age <- ifelse(pdataScoreMerge$`patient age at primary diagnosis (years):ch1`>=80,">=80",
                              ifelse(pdataScoreMerge$`patient age at primary diagnosis (years):ch1` <=40,"<=40",
                                     ifelse(pdataScoreMerge$`patient age at primary diagnosis (years):ch1` > 40 & pdataScoreMerge$`patient age at primary diagnosis (years):ch1` <=60,"40-60","60-80")))
pdataScoreMerge$Age <- factor(pdataScoreMerge$Age,levels = c("<=40","40-60","60-80",">=80"))
pdataScoreMerge$Stage <- pdataScoreMerge$`stage at primary diagnosis 5th edition:ch1`
pdataScoreMerge$Stage <- factor(pdataScoreMerge$Stage,levels = c("Stage I","Stage II","Stage III"))

cox.res=coxph(Surv(OS.time.year,OS.status)~ Group+Gender+Age+Stage,data = pdataScoreMerge)
cox.res
    
ggforest(model = cox.res,data = pdataScoreMerge,
        main = "Hazard ratio in GSE54467",
        cpositions = c(0.00, 0.20, 0.35),
        fontsize = 1.0,
        refLabel = "1", noDigits = 4)
