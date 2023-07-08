# aim: survival analysis 
# "GSE91061", "PRJEB23709", "phs000452.v3" 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/05.ICB survival and response/"))
library(dplyr)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)

ICB_CIBERSORT_PScore <- readr::read_rds("ICB_bulk_survival_response/ICB_CIBERSORT_PScore.rds.gz")
ICB_plot_OS_used <- ICB_CIBERSORT_PScore[ICB_CIBERSORT_PScore$GEO_ID %in% c("GSE91061", "PRJEB23709", "phs000452.v3"),]


# ============================================================================ #
# survival ####
# ============================================================================ #

# Unify survival time into years

# PRJEB23709: month  REF: Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy.
# GSE91061: month
# phs000452.v3: day
for(i in 1:3){
  pdata_temp <- ICB_plot_OS_used$Pyroptosis_score_pdata[[i]]
  if(i == 3){
    pdata_temp$OS_Time_use <- pdata_temp$OS_Time_use/365
  } else {
    pdata_temp$OS_Time_use <- pdata_temp$OS_Time_use/12
  }
  ICB_plot_OS_used$Pyroptosis_score_pdata[[i]] <- pdata_temp
}


prop_temp = 0.25
PScore_method = "PRGs_prot" 

ICB_plot_OS_used <- ICB_plot_OS_used %>% 
  dplyr::mutate(
    # Pre-treatment
    PScore_pdata_pre = purrr::map2(.x=Pyroptosis_score_pdata,.y = GEO_ID,function(.x,.y){
      if (length(grep("PrePost_use",colnames(.x))) > 0) {.x <- .x[.x$PrePost_use %in% "PRE",]}
      Score_pdata <- .x
      # add group information -- surv_cutpoint
      res.cut <- surv_cutpoint(Score_pdata, time = "OS_Time_use", event = "OS_use",variables = PScore_method,minprop = prop_temp)
      res.cat <- surv_categorize(res.cut)
      Score_pdata$Group <- ifelse(res.cat[[PScore_method]] == "high","High_PS","Low_PS")
      Score_pdata$Group <- factor(Score_pdata$Group,levels = c("Low_PS","High_PS"))
      return(Score_pdata)}),
    # survival analysis
    survival_plot = purrr::map2(.x=PScore_pdata_pre,.y = GEO_ID,function(.x,.y){
      Score_pdata <- .x
      
      diff=survdiff(Surv(OS_Time_use, OS_use) ~ Group,data = Score_pdata)
      pValue=1-pchisq(diff$chisq,df=1)
      pValue=signif(pValue,2)
      print(pValue)
      pValue=format(pValue, scientific = TRUE)
      
      cox.res=coxph(Surv(OS_Time_use, OS_use)~ Group,data = Score_pdata)
      HR <- signif(exp(coef(cox.res)),2)
      CI <- signif(exp(confint(cox.res)),2)
      CI_low <- signif(CI[1],2)
      CI_high <- signif(CI[2],2)
      
      Score_pdata$OS_year <- Score_pdata$OS_Time_use
      fit <- survfit(Surv(OS_year, OS_use) ~ Group, data = Score_pdata)
      
      survival_plot <- ggsurvplot(fit, 
                                  data=Score_pdata,
                                  conf.int=F,
                                  pval=paste0("P = ",pValue,"\nHR = ",paste0(HR," (",CI_low," - ",CI_high,")")),
                                  pval.size=4,
                                  risk.table=TRUE,
                                  legend.labs =  names(table(Score_pdata$Group)) ,
                                  legend.title="Cluster",
                                  title=.y,
                                  xlab="Time(years)",
                                  ylab="Overall Survival",
                                  risk.table.title="",
                                  palette=c("#0072B5FF","#BC3C29FF"),  
                                  risk.table.height=.3)
      return(survival_plot)}))


save(ICB_plot_OS_used, file = "ICB_bulk_survival_response/ICB_plot_OS_used.GSE91061_PRJEB23709_phs000452.v3.RData")

