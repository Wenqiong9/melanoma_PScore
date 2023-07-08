

rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/02.mutation/"))
library(survival)
library(survminer)

load(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/PScoreModel_ROC/PyropScore_SKCM.Rdata"))
load("D:/Work/Data/TCGAmutations/extdata/MC3/SKCM.RData")
clinicalData <- read.delim("D:/Work/Data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
clinicalData <- clinicalData[,c(1,4,6,22,25,26)]
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/heatmap/PScore_group_meta.RData")


# ============================================================================= #
# data prepare
# ============================================================================= #
patientCluster <- list("SKCM (Metastatic)" = data_ES_List$metastatic$sample_id)
gene = "BRAF"

SKCM <- tcga_skcm_mc3
SKCM_mut <- SKCM@data
SKCM_mut <-SKCM_mut[SKCM_mut$Hugo_Symbol == gene,]
MutSample <- unique(SKCM_mut[!SKCM_mut$Variant_Classification=="Nonsense_Mutation",]$Tumor_Sample_Barcode)
Sample <- SKCM@variants.per.sample$Tumor_Sample_Barcode

clinicalData$OS.time <- as.numeric(as.character(clinicalData$OS.time))
clinicalData$OS.year <- clinicalData$OS.time/365
clinicalData$OS <- as.numeric(as.character(clinicalData$OS))
clinicalData <- clinicalData[clinicalData$bcr_patient_barcode %in% Sample,]


# ============================================================================= #
# braf survival
# ============================================================================= #
clinicalData$Mut_type <- ifelse(clinicalData$bcr_patient_barcode %in% MutSample,gene,"WT")

Cluster_metastatic <- gsub("\\.","-",substr(patientCluster$`SKCM (Metastatic)`,1,12))
Clinical_metastatic <- clinicalData[clinicalData$bcr_patient_barcode %in% Cluster_metastatic,]

fit <- survfit(Surv(OS.year, OS) ~ Mut_type, data = Clinical_metastatic)

p1 <- ggsurvplot(fit,risk.table=TRUE, 
                  conf.int=F,
                  palette = "nejm",
                  pval=T,
                  pval.size=4,
                  pval.method=T,
                  legend.title="",
                  xlab="Time(years)",
                  ylab="Overall Survival",)
p1

# pdf(file="mutation_BRAF_survival_metastatic.pdf",width=6,height = 6)
# p1
# dev.off()


# ============================================================================= #
# braf x PScore_group (High_PS_BRAF,High_PS_WT,Low_PS_BRAF,Low_PS_WT)
# ============================================================================= #

data_ES <- data_ES_List$metastatic
data_ES$barcode <- substr(gsub("\\.","-",rownames(data_ES)),1,12)
melcli_group <- merge(data_ES,Clinical_metastatic,by.x="barcode",by.y="bcr_patient_barcode",all.x = TRUE,all.y = TRUE)
melcli_group <- melcli_group[which(melcli_group$OS!="" & melcli_group$OS.time!="" ),]

# PScore_group
# melcli_group$Group <- ifelse(melcli_group$PRGs_prot > unique(PScore_group_meta$cutpoint),"High_PS","Low_PS")
melcli_group <- merge(melcli_group, PScore_group_meta, by = "PRGs_prot")

# PScore_group x mutation
melcli_group$Group <- paste0(melcli_group$Group,"_",melcli_group$Mut_type)
table(melcli_group$Group)

diff=survdiff(Surv(OS.year, OS) ~ Group,data = melcli_group)
pValue=1-pchisq(diff$chisq,df=1)
pValue=format(signif(pValue,2), scientific = TRUE)
print(paste0("p value is ",pValue))

fit <- survfit(Surv(OS.year, OS) ~ Group, data = melcli_group)
p2 <- ggsurvplot(fit, 
                 conf.int=F,
                 pval=paste0("p=",pValue),
                 pval.size=4,
                 risk.table=TRUE,
                 # legend.labs=c("High_PS_BRAF", "High_PS_WT", "Low_PS_BRAF", "Low_PS_WT"),
                 legend.title="",
                 xlab="Time(years)",
                 ylab="Overall Survival",
                 palette=c("nejm"))
p2

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig2c-d.pdf",onefile = FALSE,width = 12, height = 6)
arrange_ggsurvplots(list(p1,p2),ncol=2,nrow = 1)
dev.off()



# ============================================================================= #
# braf x PScore_group forest plot
# ============================================================================= #

melcli_group$Group <- factor(melcli_group$Group,levels = c("Low_PS_WT","Low_PS_BRAF","High_PS_BRAF","High_PS_WT"))

cox.res=coxph(Surv(OS.year, OS)~ Group,data = melcli_group)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/figS5e.pdf",onefile = FALSE,width = 8, height = 4)
ggforest(cox.res, data = melcli_group,
          main="Hazard ratio in patients with different status of BRAF mutation and PScores", 
          cpositions= c(0.10, 0.22, 0.4), 
          fontsize = 1.0, refLabel = "reference", 
          noDigits= 4)
dev.off()
