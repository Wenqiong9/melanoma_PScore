# aim: response
# "GSE35640", "PRJEB23709"
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/05.ICB survival and response/"))
library(dplyr)
library(ggplot2)
library(ggpubr)

ICB_CIBERSORT_PScore <- readr::read_rds("ICB_bulk_survival_response/ICB_CIBERSORT_PScore.rds.gz")
ICB_CIBERSORT_Response <- ICB_CIBERSORT_PScore[c(1,8),]
ICB_CIBERSORT_Response$GEOID_Cancer <- paste(ICB_CIBERSORT_Response$GEO_ID,ICB_CIBERSORT_Response$cancer,sep="_")


# ============================================================================ #
# boxplot ####
# ============================================================================ #
PScore_method = "PRGs_prot"

ICB_Response_Plot <-  purrr::map2(.x=ICB_CIBERSORT_Response$Pyroptosis_score_pdata,.y=ICB_CIBERSORT_Response$GEOID_Cancer, function(.x,.y){
  
  print(paste0("Running for ",.y,"..."))
  if (length(grep("PrePost_use",colnames(.x))) > 0) {.x <- .x[.x$PrePost_use %in% "PRE",]}
  
  Score_pdata <- .x
  Score_pdata <- Score_pdata[!is.na(Score_pdata$response_use),]  # GSE35640
  Score_pdata$response_use <- factor(Score_pdata$response_use,levels = c("Benefit","NonBenefit"))
  Score_pdata$PyropScore <- scale(Score_pdata[[PScore_method]])
  
  compare <- list(c("NonBenefit","Benefit"))
  pp <- ggplot(Score_pdata,aes(x=response_use,y=PyropScore,color=response_use))+ #,fill=Benefit
    geom_boxplot(width = 0.6)+
    ggbeeswarm::geom_quasirandom(aes(color=response_use),size=2)+
    scale_x_discrete(limit=c("NonBenefit","Benefit"),
                     labels=c(paste0("NR (n=",table(Score_pdata$response_use)[2],")"),
                              paste0("R (n=",table(Score_pdata$response_use)[1],")")))+
    scale_fill_manual(values = c("#E41A1C","#377EB8"),name="Response")+
    scale_color_manual(values = c("#E41A1C","#377EB8"),guide = "none")+
    stat_compare_means(comparisons = compare)+
    xlab("") + ylab("Pyroptosis Score")+ labs(title = .y)+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_text(color="black",size=12,angle = 60,hjust = 1,vjust = 1),  #angle = 60 横坐标角度
            axis.text.y = element_text(color="black",size=10),
            axis.title.y = element_text(color="black",size=12),
            axis.ticks.x = element_blank(),legend.position = "top",legend.direction = "horizontal",
            axis.line.y = element_line(colour = "black"),
            axis.line.x = element_line(colour = "black"),
            strip.background = element_blank())
  return(pp)
})

# pdf("ICB_bulk_survival_response/ICB_bulk_response_PRCR_PDSD.pdf", width = 20, height = 8)
cowplot::plot_grid(plotlist = ICB_Response_Plot,ncol=2)
# dev.off()
  

# ============================================================================ #
# chisq test ####
# ============================================================================ #
PScore_method = "PRGs_prot"  

ICB_chisq_Plot <-  purrr::map2(.x=ICB_CIBERSORT_Response$Pyroptosis_score_pdata,.y=ICB_CIBERSORT_Response$GEO_ID, function(.x,.y){
  
  print(paste0("Running for ",.y,"..."))
  if (length(grep("PrePost_use",colnames(.x))) > 0) {.x <- .x[.x$PrePost_use %in% "PRE",]}
  immune_dataset <- .x
  immune_dataset <- immune_dataset[immune_dataset$response_use!="None",]
  immune_dataset$Benefit <- immune_dataset$response_use 
  
  # group
  if(.y == "GSE35640"){
    immune_dataset$Group2 <- ifelse(immune_dataset$PRGs_prot >= quantile(immune_dataset$PRGs_prot,0.55), "High_PS", "Low_PS")
  } else {
    res.cut <- surv_cutpoint(immune_dataset, time = "OS_Time_use", event = "OS_use",variables = "PRGs_prot", minprop = 0.2)
    res.cat <- surv_categorize(res.cut)
    immune_dataset$Group2 <- res.cat[[PScore_method]]
    immune_dataset$Group2 <- ifelse(immune_dataset$Group2=="high","High_PS","Low_PS")
  }
  immune_dataset$Group2 <- factor(immune_dataset$Group2, levels = c("Low_PS","High_PS"))
  
  tmp <- table(immune_dataset$Group2,immune_dataset$Benefit)
  if(min(tmp)>=5){
    chisq_p <- chisq.test(immune_dataset$Group2,immune_dataset$Benefit)
    pvalue <- chisq_p$p.value
  }else{
    fisher_p <- fisher.test(immune_dataset$Group2,immune_dataset$Benefit)
    pvalue <- fisher_p$p.value
  }
  
  pp <- ggplot(data=immune_dataset, mapping=aes(x=Group2,fill=Benefit))+
    geom_bar(stat="count",width=0.5,position='fill')+
    scale_fill_manual(values=c("#BC3C29FF","#0072B5FF"))+  
    geom_text(stat='count',aes(label=..count..) 
              , color="white", size=3.5,position=position_fill(0.5))+
    theme(panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(color="black",size=12,angle = 60,hjust = 1,vjust = 1),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "top",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          strip.background = element_blank())+
    labs(title=.y, subtitle = paste0("p value = ",format(signif(pvalue,2),scientific = T)) , y="Pecentage(%)")
  pp
  return(pp)})

cowplot::plot_grid(plotlist = ICB_chisq_Plot,ncol=2)

save(ICB_Response_Plot,ICB_chisq_Plot,file="d:/Work/Projects/BMC_Revised/Results/06.figures/fig6b-c_PRJEB23709_GSE35640.RData")
