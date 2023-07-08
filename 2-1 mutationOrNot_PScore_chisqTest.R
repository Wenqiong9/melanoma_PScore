# braf chisq test


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/02.mutation/"))

load("D:/Work/Data/TCGAmutations/extdata/MC3/SKCM.RData")
SKCM <- tcga_skcm_mc3@data
SKCM_mut<-SKCM[SKCM$Hugo_Symbol == "BRAF",]

load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/heatmap/PScore_group_meta.RData") 
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/PScoreModel_ROC/PyropScore_SKCM.Rdata") 


# ============================================================================ #
# data prepare
# ============================================================================ #
data_ES_m <- data_ES_List$metastatic


# merge BRAF status
data_ES_m$barcode <- substr(gsub("\\.", "-", rownames(data_ES_m)),1,16)
mut_total <- substr(unique(SKCM$Tumor_Sample_Barcode_full),1,16)
mut_BRAF <- substr(unique(SKCM_mut$Tumor_Sample_Barcode_full),1,16)
# mut_BRAF <- substr(unique(SKCM_mut[!SKCM_mut$Variant_Classification=="Nonsense_Mutation",]$Tumor_Sample_Barcode_full),1,16)
## patient TCGA-D9-A6EC who had Nonsense_Mutation also had Missense_Mutation
data_ES_m$BRAF_status <- ifelse(data_ES_m$barcode %in% mut_BRAF, "BRAF", ifelse(data_ES_m$barcode %in% mut_total, "WT", NA))
save(data_ES_m,file = "mutation_complexHeatmap.RData")

# check data
table(data_ES_m$barcode %in% mut_total)
# FALSE  TRUE 
# 5   364 
table(data_ES_m$barcode %in% mut_BRAF)
# FALSE  TRUE 
# 180   189
table(data_ES_m$BRAF_status)
# BRAF   WT 
# 189  175 


# merge PScore group
data_ES_m <- merge(data_ES_m, PScore_group_meta, by = "PRGs_prot")
# data_ES_m$group <- ifelse(data_ES_m$PRGs_prot > unique(PScore_group_meta$cutpoint), "High_PS", "Low_PS")


# ============================================================================ #
# chisq_p
# ============================================================================ #
library(ggplot2)

chisq_p <- chisq.test(data_ES_m$Group,data_ES_m$BRAF_status)
pvalue <- format(signif(chisq_p$p.value,2), scientific = TRUE)

p <- 
ggplot(data=data_ES_m, mapping=aes(x=Group,fill=BRAF_status))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c("#A73030FF","#A7303033"))+ 
  geom_text(stat='count',aes(label= after_stat(count)) 
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
  labs(title = paste0("p value = ",pvalue), y="Pecentage(%)")
p

library(patchwork)
load("braf_boxplot.RData")

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig2b.pdf", width = 10, height = 7)
braf_boxplot|p
dev.off()

