# TCGA surv_cutpoint survival analysis & multiCOX
# 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
library(survival)
library(survminer)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))

load("PScoreModel_ROC/PyropScore_SKCM.Rdata")
clinic <- read.delim("d:/Work/Data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
mel_cli <- clinic[clinic$type%in%"SKCM",]  


# ============================================================================ #
# 1. 生存分析####
# ============================================================================ #
PScore_method = names(data_ES_List$metastatic)[3] # PRGs_prot
 
# merge PScore and clinical data
data_ES <- data_ES_List$metastatic
data_ES$barcode <- rownames(data_ES)
data_ES$barcode <- gsub("\\.","_",substr(data_ES$barcode,1,12))
mel_cli$bcr_patient_barcode <- gsub("-","_",mel_cli$bcr_patient_barcode)
melcli_group <- merge(data_ES,mel_cli,by.x="barcode",by.y="bcr_patient_barcode",all.x = TRUE)

# survival data prepare
melcli_group$OS <- as.numeric(melcli_group$OS)
melcli_group$OS.time <- as.numeric(melcli_group$OS.time)
melcli_group <- melcli_group[which(melcli_group$OS!="" & melcli_group$OS.time!="" ),] 
melcli_group$OS_year <- melcli_group$OS.time/365

# cutpoint 
res.cut <- surv_cutpoint(melcli_group, time = "OS.time", event = "OS",variables = PScore_method, minprop = 0.4)
res.cat <- surv_categorize(res.cut)
melcli_group$Group <- ifelse(res.cat[[PScore_method]] == "high","High_PS","Low_PS")
melcli_group$Group <- factor(melcli_group$Group,levels = c("Low_PS","High_PS"))

diff=survdiff(Surv(OS_year, OS) ~ Group,data = melcli_group)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue=format(pValue, scientific = TRUE)
pValue

# cox HR
cox.res=coxph(Surv(OS_year, OS)~ Group,data = melcli_group)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)
CI_low <- signif(CI[1],2)
CI_high <- signif(CI[2],2)

fit <- survfit(Surv(OS_year, OS) ~ Group, data = melcli_group)
    
ggsurvplot(fit, 
           data=melcli_group,
           conf.int=F,
           pval=paste0("p = ",pValue,"\nHR = ",paste0(HR,"(",CI_low," - ",CI_high,")")),
           risk.table=TRUE,
           legend.labs=c("Low_PS","High_PS"),
           legend.title="") 
          

# ============================================================================ #
# 2.multicox####
# ============================================================================ #

# clinical factors
melcli_group$age <- ifelse(melcli_group$age_at_initial_pathologic_diagnosis>=80,">=80",
                           ifelse(melcli_group$age_at_initial_pathologic_diagnosis<=40,"<=40",
                                  ifelse(melcli_group$age_at_initial_pathologic_diagnosis>40&melcli_group$age_at_initial_pathologic_diagnosis<=60,"40-60","60-80")))
melcli_group$age <- factor(melcli_group$age,levels = c("<=40","40-60","60-80",">=80"))

melcli_group$stage <- melcli_group$ajcc_pathologic_tumor_stage
melcli_group$stage <- gsub("A|B|C|D","",melcli_group$stage)
melcli_group$stage <- factor(melcli_group$stage,levels = c("Stage I","Stage II","Stage III","Stage IV"))

melcli_group$gender[melcli_group$gender=="MALE"] <- "Male"
melcli_group$gender[melcli_group$gender=="FEMALE"] <- "Female"
melcli_group$gender <- factor(melcli_group$gender,levels = c("Male","Female"))


# forest plot
library(survival)
library(survminer)
cox.res=coxph(Surv(OS_year,OS)~ Group+gender+age+stage,data = melcli_group)
cox.res

ggforest(
  model = cox.res,data = melcli_group,
  main = "Hazard ratio in TCGA-SKCM",
  cpositions = c(0.00, 0.20, 0.35),
  fontsize = 1.0,
  refLabel = "1", noDigits = 4)

PScore_group_meta <- melcli_group[,c("barcode", "PRGs_prot", "Group")]
PScore_group_meta$cutpoint <- res.cut$cutpoint$cutpoint
save(PScore_group_meta,file = "heatmap/PScore_group_meta.RData")
