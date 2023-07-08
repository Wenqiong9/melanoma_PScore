# inhouse data
# 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/05.ICB survival and response/"))
library(dplyr)
library(ggplot2)
load((paste0(path, "Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")))
inhouse_expr <- openxlsx::read.xlsx("d:/Work/Data/inhouse_data(xiangya)/Supplementary tables_HEYi_Ferroptosis.xlsx",sheet=9,startRow = 2,rowNames = TRUE)
inhouse_meta <- openxlsx::read.xlsx("d:/Work/Data/inhouse_data(xiangya)/Supplementary tables_HEYi_Ferroptosis.xlsx",sheet=8,startRow = 2)

genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])


# ============================================================================ #
# 1.data prepare -- inhouse data####
# ============================================================================ #

# ssgsea
genelist$PRGs_prot[genelist$PRGs_prot %in% rownames(inhouse_expr)==FALSE]
rownames(inhouse_expr)[rownames(inhouse_expr) == "PJVK"] <- "DFNB59"

expr_Data <- inhouse_expr %>% as.matrix()
data_ES <- GSVA::gsva(expr_Data,genelist,method="ssgsea")
data_ES <- data_ES %>% t() %>% data.frame()
head(data_ES)

# merge PScore and response information
data_ES$ID <- rownames(data_ES)
inhouse_meta <- merge(data_ES,inhouse_meta,by="ID")
save(inhouse_meta, file = "ICB_bulk_survival_response/inhouse_meta_PScore.RData")


# ============================================================================ #
# 2.boxplot show the relationship between PScore and response ####
# ============================================================================ #
PScore_method = "PRGs_prot"
compare <- list(c("R","NR"))

boxplot_inhouse <- ggplot(inhouse_meta,aes(x=Response,y=.data[[PScore_method]],color=Response))+ #,fill=Benefit
    geom_boxplot(width = 0.6)+
    ggbeeswarm::geom_quasirandom(aes(color=Response),size=2)+
    scale_x_discrete(limit=c("NR","R"),
                     labels=c(paste0("NR (n=",sum(inhouse_meta$Response=="NR"),")"),
                              paste0("R (n=",sum(inhouse_meta$Response=="R"),")")))+
    scale_fill_manual(values = c("#377EB8","#E41A1C"),name="Response")+
    scale_color_manual(values = c("#377EB8","#E41A1C"),guide = "none")+
    ggpubr::stat_compare_means(comparisons = compare)+
    xlab("") + ylab("Pyroptosis Score")+ labs(title = paste0("inhouse data ",PScore_method))+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_text(color="black",size=12,angle = 60,hjust = 1,vjust = 1),  
            axis.text.y = element_text(color="black",size=10),
            axis.title.y = element_text(color="black",size=12),
            axis.ticks.x = element_blank(),legend.position = "top",legend.direction = "horizontal",
            axis.line.y = element_line(colour = "black"),
            axis.line.x = element_line(colour = "black"),
            strip.background = element_blank())

load("d:/Work/Projects/BMC_Revised/Results/06.figures/fig6b-c_PRJEB23709_GSE35640.RData")
pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig6b.pdf", width = 12,height = 6)
cowplot::plot_grid(plotlist = list(boxplot_inhouse, ICB_Response_Plot[[1]], ICB_Response_Plot[[2]]),ncol = 3)
dev.off()

# ============================================================================ #
# 3.survival analysis####
# ============================================================================ #
library(survminer)
library(survival)
inhouse_meta$Response <- factor(inhouse_meta$Response,levels = c('R','NR'))

PScore_method = "PRGs_prot"  
prop_temp = 0.2

res.cut <- surv_cutpoint(inhouse_meta, time = "PFST.Month.", event = "PFS",variables = "PRGs_prot", minprop = prop_temp)
res.cat <- surv_categorize(res.cut)
inhouse_meta$Group <- ifelse(res.cat[[PScore_method]] == "high","High_PS","Low_PS")
inhouse_meta$Group <- factor(inhouse_meta$Group,levels = c("Low_PS","High_PS"))

inhouse_meta$PFS <- as.numeric(inhouse_meta$PFS )
inhouse_meta$PFST.Month. <- as.numeric(inhouse_meta$PFST.Month. )

diff=survdiff(Surv(PFST.Month., PFS) ~ Group,data = inhouse_meta)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue
pValue=format(pValue, scientific = TRUE)

cox.res=coxph(Surv(PFST.Month., PFS)~ Group,data = inhouse_meta)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)
CI_low <- signif(CI[1],2)
CI_high <- signif(CI[2],2)

inhouse_meta$PFS_year <- inhouse_meta$PFST.Month./12
fit <- survfit(Surv(PFS_year, PFS) ~ Group, data = inhouse_meta)

survival_plot <- ggsurvplot(fit, 
                            data=inhouse_meta,
                            conf.int=F,
                            pval=paste0("P = ",pValue,"\nHR = ",paste0(HR," (",CI_low," - ",CI_high,")")),
                            pval.size=4,
                            risk.table=TRUE,
                            legend.labs =  names(table(inhouse_meta$Group)) ,
                            legend.title="Cluster",
                            title="inhouse_data",
                            xlab="Time(years)",
                            ylab="PFS",
                            risk.table.title="",
                            palette=c("#0072B5FF","#BC3C29FF"), 
                            risk.table.height=.3)

# add TCGA and other datasets -- 5-2.R ##
load("ICB_bulk_survival_response/ICB_plot_OS_used.GSE91061_PRJEB23709_phs000452.v3.RData")  # ICB_plot_OS_used
load("ICB_bulk_survival_response/ICB_plot_OS_used.TCGA.RData")  # p

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig6d.pdf",width=18,height=5,onefile = FALSE)
arrange_ggsurvplots(list(survival_plot, ICB_plot_OS_used$survival_plot[[1]], ICB_plot_OS_used$survival_plot[[2]],
                         ICB_plot_OS_used$survival_plot[[3]], p),ncol=5)
dev.off()


# ============================================================================ #
# 4.chisq test ####
# ============================================================================ #
fisher_p <- fisher.test(inhouse_meta$Group,inhouse_meta$Response)
fisher_p <- format(signif(fisher_p$p.value,2), scientific = TRUE)

p_barplot <- 
  ggplot(data=inhouse_meta, mapping=aes(x=Group,fill=Response))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c("#BC3C29FF","#0072B5FF"))+  #"#377EB8","#E41A1C"  金色和灰色"#FFC000","#A5A5A5"
  geom_text(stat='count',aes(label=..count..) #scales::percent(..count../sum(..count..)
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
  labs(title="inhouse_data", subtitle = paste0("p value = ",fisher_p) , y="Pecentage(%)")


pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig6c.pdf", width = 12,height = 6)
cowplot::plot_grid(plotlist = list(p_barplot, ICB_chisq_Plot[[1]], ICB_chisq_Plot[[2]]),ncol = 3)
dev.off()


cowplot::plot_grid(plotlist = list(p_barplot), ncol=1, nrow = 1)
dev.off()


# ============================================================================ #
# 5.inhouse data heatmap ####
# ============================================================================ #

library(ComplexHeatmap)
# load("ICB_bulk_survival_response/inhouse_meta_PScore.RData")
# inhouse_expr <- openxlsx::read.xlsx("d:/Work/Data/inhouse_data(xiangya)/Supplementary tables_HEYi_Ferroptosis.xlsx",sheet=9,startRow = 2,rowNames = TRUE)

# main of heatmap
melanomaExpr <- inhouse_expr[rownames(inhouse_expr) %in% genelist$PRGs_prot,]
melanomaExpr <- t(scale(t(melanomaExpr)))

# patient order of patients
inhouse_meta <- inhouse_meta[order(inhouse_meta[[PScore_method]]),]
melanomaExpr <- melanomaExpr[,match(inhouse_meta$ID,colnames(melanomaExpr))]

# anno
heatmapAnno <- HeatmapAnnotation(
  PFS = anno_barplot(inhouse_meta$PFST.Month.),
  Age = inhouse_meta$Age,
  Sex = inhouse_meta$Sex,
  SubTypes = inhouse_meta$SubTypes,
  Response = inhouse_meta$Response,
  Group = inhouse_meta$Group,
  PScore = inhouse_meta[[PScore_method]],
  col = list(Age = circlize::colorRamp2(c(11, 100), c("white", "#1f78b4")),
             Sex = c("Male" = "#E18727FF", "Female" = "#FFDC91FF"),
             SubTypes = c("Acral Melanoma" = "#BC3C29FF", "Cutaneous Melanoma" = "#0072B5FF",
                          "Mucosal Melanoma" = "#E18727FF", "Nodular Melanoma" = "#FFDC91FF",
                          "Unknow" = "#20854EFF"),
             Response = c("R" = "darkgreen", "NR" = "lightgreen"),
             Group = c("High_PS"="#BC3C29FF", "Low_PS"="#0072B5FF"),
             PScore =  circlize::colorRamp2(c(1.5, 2.5), c("white", "#A73030FF"))))

# heatmapAnnoRow <- rowAnnotation(
#   Prognosis = ifelse(rownames(melanomaExpr) %in% uniCoxSig_list$metastatic_OS$gene[uniCoxSig_list$metastatic_OS$Risk_group == "Protect"], "Protect",
#                      ifelse(rownames(melanomaExpr) %in% uniCoxSig_list$metastatic_OS$gene[uniCoxSig_list$metastatic_OS$Risk_group == "Risk"],"Risk","Not significant")),
#   col = list(Prognosis = c("Protect" = "#0072B5EE", "Risk" = "#BC3C29EE", "Not significant" = "#0072B599")))

p_heatmap <- Heatmap(melanomaExpr, 
                     rect_gp = gpar(col= "white"),
                     cluster_columns = FALSE,
                     show_column_names = FALSE,
                     show_row_dend = T, 
                     row_names_side = "left", 
                     row_km = 2, 
                     column_title="Expression of Interested Genes",
                     name =  "Expressions",
                     top_annotation = heatmapAnno) #,
                     # left_annotation = heatmapAnnoRow)
p_heatmap


pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig6a.pdf",width = 8,height = 8)
p_heatmap
dev.off()
