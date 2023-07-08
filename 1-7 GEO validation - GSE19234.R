# GEO median survival analysis & multiCOX
# GSE19234 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggplot2)
library(ggpubr)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
load("coxForest/uniCoxSig_list.RData")
load("d:/Work/Data/GEO/GEO_melanoma_res/Melanoma_Datasets_16.RData")

genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])


# ============================================================================ #
# 1.data prepare ####
# ============================================================================ #
ExprMatrix <- Melanoma_Datasets_16$gset_data_matrix_symbol[[4]]
ExprMatrix <- log2(ExprMatrix+1) 

genelist$PRGs_prot[which(genelist$PRGs_prot %in% rownames(ExprMatrix)==FALSE)]
rownames(ExprMatrix)[rownames(ExprMatrix)=="PJVK"] <- "DFNB59"


pData <- Melanoma_Datasets_16$gset_data_pData[[4]]
pData <- pData[,c(1:2,(ncol(pData)-5):ncol(pData))]

# delete duplicated samples
pData <- pData[-which(nchar(strsplit(pData$title,split = ""))==16),]

# delete stage IV, because multivariate cox delete this 
pData <- pData[pData$`stage at reccurence:ch1` != 'IV',]
# table(pData$`stage at reccurence:ch1`)
# IIIA IIIB IIIC   IV 
#    4   18   11    5 

# ============================================================================ #
# 2.PScore ####
# ============================================================================ #
data_ES <- GSVA::gsva(as.matrix(ExprMatrix),genelist,method="ssgsea")
data_ES <- as.data.frame(t(data_ES))
data_ES_zscore <- scale(data_ES) %>% data.frame()  


load("d:/Work/Projects/BMC_Revised/Results/07.tables/PScores_of_patients.RData")
data_ES_save <- data_ES
data_ES_save$patient_id <- rownames(data_ES_save)
data_ES_save$dataset <- "GSE19234"
PScores_of_patients[["GSE19234"]] <- data_ES_save
save(PScores_of_patients, file = "d:/Work/Projects/BMC_Revised/Results/07.tables/PScores_of_patients.RData")


data_ES$ID <- rownames(data_ES)
gse_pdata_m <- merge(data_ES,pData,by.x="ID",by.y="geo_accession")


# ============================================================================ #
# 3.survival ####
# ============================================================================ #
library(survminer)
library(survival)

gse_pdata_m$OS.status <- gse_pdata_m$`staus dead or alive:ch1` %>% as.numeric()
gse_pdata_m$OS.time <- gse_pdata_m$`days since initial diagnosis:ch1` %>% as.numeric()

# group
gse_pdata_m$PScore <- gse_pdata_m$PRGs_prot
res.cut <- surv_cutpoint(gse_pdata_m, time = "OS.time", event = "OS.status",variables = "PScore",minprop = 0.4)
res.cat <- surv_categorize(res.cut)
gse_pdata_m$Group <- ifelse(res.cat$PScore == "high","High_PS","Low_PS")
gse_pdata_m$Group <- factor(gse_pdata_m$Group,levels = c("Low_PS","High_PS"))
table(gse_pdata_m$Group)

# p 
diff=survdiff(Surv(OS.time,OS.status) ~ Group,data = gse_pdata_m)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue=format(pValue, scientific = TRUE)
pValue

# cox
cox.res=coxph(Surv(OS.time,OS.status)~ Group,data = gse_pdata_m)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)
CI_low <- signif(CI[1],2)
CI_high <- signif(CI[2],2)

gse_pdata_m$OS.time.year <- gse_pdata_m$OS.time/365
fit <- survfit(Surv(OS.time.year,OS.status) ~ Group, data = gse_pdata_m)

ggsurvplot(fit, 
           data=gse_pdata_m,
           conf.int=F,
           pval=paste0("p=",pValue,"\nHR=",paste0(HR,"(",CI_low,"-",CI_high,")")),
           risk.table=TRUE,
           legend.labs=factor(c("Low_PS","High_PS")),
           legend.title="")


# ================ #
# 4.forest plot ####
# ================ #

# clinical factors
gse_pdata_m$Gender <- gse_pdata_m$`gender:ch1`
gse_pdata_m$Gender[gse_pdata_m$Gender=="male"] <- "Male"
gse_pdata_m$Gender[gse_pdata_m$Gender=="female"] <- "Female"
gse_pdata_m$Gender <- factor(gse_pdata_m$Gender,levels = c("Male","Female"))

gse_pdata_m$Age <- gse_pdata_m$`age:ch1` %>% as.numeric()
gse_pdata_m$Age <- ifelse(gse_pdata_m$Age>=80,">=80",ifelse(gse_pdata_m$Age<=40,"<=40",ifelse(gse_pdata_m$Age>40&gse_pdata_m$Age<=60,"40-60","60-80")))
gse_pdata_m$Age <- factor(gse_pdata_m$Age,levels = c("<=40","40-60","60-80",">=80"))

gse_pdata_m$Stage <-  gse_pdata_m$`stage at reccurence:ch1`
gse_pdata_m$Stage <- factor(gse_pdata_m$Stage)

cox.res=coxph(Surv(OS.time.year,OS.status)~ Group+Gender+Age,data = gse_pdata_m)  # +Stage 
cox.res
# caused by too many variables, need to re-filter your variables
# Warning message:
#   In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  :
#                  Loglik converged before variable  6,7 ; coefficient may be infinite.
    
ggforest(model = cox.res,data = gse_pdata_m,
          main = "Hazard ratio in GSE19234",
          cpositions = c(0.00, 0.20, 0.35),
          fontsize = 1.0,
          refLabel = "1", noDigits = 4)

