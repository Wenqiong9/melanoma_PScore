# Aim: survival of TCGA-SKCM-ICB 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/05.ICB survival and response/"))
library(survival)
library(survminer)
load("d:/Work/Data/TCGA_secondary/data_clinical_clean/data_clinical_merge.RData")
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/PScoreModel_ROC/PyropScore_SKCM.Rdata")


# ============================================================================ # 
# 1.data prepare ####
# ============================================================================ #
# check data
dim(data_clinical_merge)
# [1] 523  20
table(data_clinical_merge$bcr_patient_barcode %in% substr(rownames(data_ES_List$metastatic),1,12))
table(data_clinical_merge$bcr_patient_barcode %in% substr(rownames(data_ES_List$primary),1,12))
# FALSE  TRUE 
# 111   412 
# FALSE  TRUE 
# 410   113 


# TCGA_ICB
data_clinical_merge <- data_clinical_merge[data_clinical_merge$Therapy_Type == "Immunotherapy",]
ICB_with_response <- unique(data_clinical_merge$bcr_patient_barcode[!is.na(data_clinical_merge$Drug_Response)])
table(ICB_with_response %in% substr(rownames(data_ES_List$metastatic),1,12))
table(ICB_with_response %in% substr(rownames(data_ES_List$primary),1,12))
# FALSE  TRUE 
# 4    39 
# FALSE  TRUE 
# 39     4 


# merge PScore and clinic information
data_ES_List$metastatic$bcr_patient_barcode <- substr(rownames(data_ES_List$metastatic),1,12)
ClinicalPScore_m <- merge(data_clinical_merge[,c(1:3,6:7,12,19:20)],data_ES_List$metastatic,by="bcr_patient_barcode")

dim(unique(ClinicalPScore_m[ClinicalPScore_m$Therapy_Type=="Immunotherapy",]))
dim(ClinicalPScore_m[ClinicalPScore_m$Therapy_Type=="Immunotherapy"&!is.na(ClinicalPScore_m$Drug_Response),])
# [1] 70 14
# [1] 40 14


# ============================================================================ #
#  2.survival analysis ####
# ============================================================================ #

PScore_method = "PRGs_prot"  #  PRGs_74   PRGs_34 PRGs_prot PRGs_risk PRGs_34_minus


res.cut <- surv_cutpoint(ClinicalPScore_m, time = "OS.time", event = "OS", variables = PScore_method, minprop = 0.2)
res.cat <- surv_categorize(res.cut)
ClinicalPScore_m$Group <- ifelse(res.cat[[PScore_method]] == "high","High_PS","Low_PS")
ClinicalPScore_m$Group <- factor(ClinicalPScore_m$Group,levels = c("Low_PS","High_PS"))

diff=survdiff(Surv(OS.time, OS) ~ Group,data = ClinicalPScore_m)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue
pValue=format(pValue, scientific = TRUE)

cox.res=coxph(Surv(OS.time,OS)~ Group,data = ClinicalPScore_m)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)
CI_low <- signif(CI[1],2)
CI_high <- signif(CI[2],2)

ClinicalPScore_m$OS_year <- ClinicalPScore_m$OS.time/365
fit <- survfit(Surv(OS_year, OS) ~ Group, data = ClinicalPScore_m)

p <- ggsurvplot(fit, 
               data=ClinicalPScore_m,
               conf.int=F,
               pval=paste0("P = ",pValue,"\nHR = ",paste0(HR," (",CI_low," - ",CI_high,")")),
               pval.size=4,
               risk.table=TRUE,
               legend.labs=c("Low_PS","High_PS"),
               title = "TCGA",
               legend.title="",
               xlab="Time(years)",
               ylab="Overall Survival",
               break.time.by = 5,
               risk.table.title="",
               palette=c("#0072B5FF","#BC3C29FF"), 
               risk.table.height=.3)
p
save(p,file = "ICB_bulk_survival_response/ICB_plot_OS_used.TCGA.RData")

