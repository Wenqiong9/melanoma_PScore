

rm(list=ls())
options(stringsAsFactors = F)
library(survival)
library(dplyr)
library(GSVA)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
load("coxForest/uniCoxSig_list.RData")
load("d:/Work/Data/GEO/GEO_melanoma_res/Melanoma_Datasets_16.RData")
genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])


# ============================================================================ #
# prepare ####
# ============================================================================ #
pData <- Melanoma_Datasets_16$gset_data_pData[[13]]

ExprMatrix <- Melanoma_Datasets_16$gset_data_matrix_symbol[[13]]
ExprMatrix <- log2(ExprMatrix+1) 
genelist$PRGs_prot[which(genelist$PRGs_prot %in% rownames(ExprMatrix)==FALSE)]
rownames(ExprMatrix)[ rownames(ExprMatrix)=="PJVK"] <- "DFNB59"

data_ES <- GSVA::gsva(as.matrix(ExprMatrix),genelist,method="ssgsea")
data_ES <- t(data_ES) %>% data.frame()
data_ES_zscore <- scale(data_ES) %>% data.frame()  

data_ES$ID <- rownames(data_ES)
GSE65904_m <- merge(data_ES,pData,by.x="ID",by.y="geo_accession")
GSE65904_m <- GSE65904_m[GSE65904_m$`tumor stage:ch1` != "Primary",]


# ============================================================================ #
# survival ####
# ============================================================================ #
library(survminer)
library(survival)
library(ggplot2)

# DSS
GSE65904_m$DSS <- GSE65904_m$`disease specific survival (1=death, 0=alive):ch1` %>% as.numeric()
GSE65904_m$DSS.time <- GSE65904_m$`disease specific survival in days:ch1` %>% as.numeric()
pdataScoreMerge <- GSE65904_m
pdataScoreMerge$DSS.time.year <- pdataScoreMerge$DSS.time/365

# cutpoint group
res.cut <- surv_cutpoint(pdataScoreMerge, time = "DSS.time", event = "DSS",variables = "PRGs_prot",minprop = 0.4)
res.cat <- surv_categorize(res.cut)
pdataScoreMerge$Group <- ifelse(res.cat$PRGs_prot == "high","High_PS","Low_PS")
pdataScoreMerge$Group <- factor(pdataScoreMerge$Group,levels = c("Low_PS","High_PS"))

diff=survdiff(Surv(DSS.time,DSS) ~ Group,data = pdataScoreMerge)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue=format(pValue, scientific = TRUE)
pValue

cox.res=coxph(Surv(DSS.time,DSS)~ Group,data = pdataScoreMerge)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)
CI_low <- signif(CI[1],2)
CI_high <- signif(CI[2],2)

fit <- survfit(Surv(DSS.time.year,DSS) ~ Group, data = pdataScoreMerge)

ggsurvplot(fit, 
          data=pdataScoreMerge,
          conf.int=F,
          pval=paste0("p=",pValue,"\nHR=",paste0(HR,"(",CI_low,"-",CI_high,")")),
          risk.table=TRUE,
          legend.labs=factor(c("Low_PS","High_PS")),
          legend.title="")


# ============================================================================ #
# GSE65904: Multivariate cox regression ####
# ============================================================================ #
pdataScoreMerge$Gender <-  pdataScoreMerge$`gender:ch1`
pdataScoreMerge$Gender <- factor( pdataScoreMerge$Gender,levels = c("Male","Female"))

pdataScoreMerge$Age <-  pdataScoreMerge$`age:ch1` %>% as.numeric()
pdataScoreMerge$Age <- ifelse( pdataScoreMerge$Age>=80,">=80",ifelse( pdataScoreMerge$Age<=40,"<=40",ifelse( pdataScoreMerge$Age>40& pdataScoreMerge$Age<=60,"40-60","60-80")))
pdataScoreMerge$Age <- factor( pdataScoreMerge$Age,levels = c("<=40","40-60","60-80",">=80"))

pdataScoreMerge$Stage <-   pdataScoreMerge$`tumor stage:ch1`
pdataScoreMerge$Stage <- factor( pdataScoreMerge$Stage,levels=c("In-transit","Regional","General","Local"))
# Stage: 
# Local recurrence, In-transit metastasis(II), Regional metastasis(III), Distant metastasis(IV)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4494939/
   
# multivariate COX
cox.res=coxph(Surv(DSS.time.year,DSS)~ Group+Gender+Age+Stage,data = pdataScoreMerge)
cox.res
cox.zph(cox.res)  # p > 0.05 means we can assume that the Cox model conforms to the proportional risk assumption
#         chisq df    p
# Group  0.0126  1 0.91
# Gender 1.8001  1 0.18
# Age    0.6956  3 0.87
# Stage  3.8652  3 0.28
# GLOBAL 6.2610  8 0.62

ggforest(model = cox.res,data = pdataScoreMerge,
        main = "Hazard ratio in GSE65904", 
        cpositions = c(0.00, 0.20, 0.35),
        fontsize = 1.0,
        refLabel = "1", noDigits = 4)

