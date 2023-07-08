# boxplot
# figure S5c-d
# ============================================================================= #


rm(list=ls())
options(stringsAsFactors = F)
library(ggplot2)
library(ggpubr)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/02.mutation/"))
load("D:/Work/Data/TCGAmutations/extdata/MC3/SKCM.RData")
SKCM <- tcga_skcm_mc3
load(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/PScoreModel_ROC/PyropScore_SKCM.Rdata"))


# ============================================================================= #
# get top10 mutated genes ####
# ============================================================================= #
{geneFreqDf <- data.frame() 
for(i in unique(SKCM@data$Hugo_Symbol)){
  tmp_gene <- SKCM@data[SKCM@data$Hugo_Symbol==i,]
  tmp_df <- data.frame("gene"=i,"freq"= length(unique(tmp_gene$Tumor_Sample_Barcode_full))) 
  geneFreqDf <- rbind(geneFreqDf,tmp_df)
}
geneFreqDf <- geneFreqDf[order(geneFreqDf$freq,decreasing=T),]
geneFreqDf[1:10,]}
# save(geneFreqDf,file = "geneFreqDf.Rdata")
# write.csv(geneFreqDf, file = "geneFreqDf.csv", row.names = F)
#         gene freq
# 80       TTN  371
# 603    MUC16  339
# 98     DNAH5  258
# 117     BRAF  245
# 1860    PCLO  225
# 859    LRP1B  212
# 1102   GPR98  207
# 1263 PKHD1L1  190
# 1692   DNAH7  187
# 1285   CSMD1  183


# ============================================================================= #
# boxplot ####
# ============================================================================= #
load("geneFreqDf.Rdata")

PScore_method <- "PRGs_prot"

data_ES_zscore <- data_ES_zscore_List$metastatic
data_ES_zscore$patientID <- gsub("\\.","_",substr(rownames(data_ES_zscore),1,16))

# other interseted genes
# PTEN + KIT + NARS + TERT + NLRP1/3/7 + NLRC4
c("PTEN", "KIT", "NARS", "TERT", "NLRP1","NLRP3", "NLRP7","NLRC4") %in% geneFreqDf$gene
gene_suppliment <- c("PTEN", "KIT", "NARS", "TERT", "NLRP1","NLRP3", "NLRP7","TLR4")


mutatedGeneBoxplot <- list() 
for(i in c(geneFreqDf$gene[1:10],gene_suppliment)){

  tmp <- SKCM@data[SKCM@data$Hugo_Symbol == i & SKCM@data$sample_type_description == "Metastatic",]
  tmp$patientID <- substr(gsub("-","_",tmp$Tumor_Sample_Barcode_full), 1, 16)
  data_ES_zscore$mutate <- ifelse(data_ES_zscore$patientID %in% tmp$patientID, i, "WT")
  
  if(i == "BRAF"){print(table( data_ES_zscore$mutate ))}
  
  mutatedGeneBoxplot[[i]] <- ggplot(data_ES_zscore, aes(x=mutate, y=.data[[PScore_method]], color=mutate))+
    geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+ 
    ggbeeswarm::geom_beeswarm(alpha=0.5,size=0.5)+
    labs(y="Pyroptosis_Score")+
    ggsci::scale_color_nejm(limits=names(table(data_ES_zscore$mutate)),label=paste0(names(table(data_ES_zscore$mutate))),name="")+
    ggpubr::stat_compare_means() + 
    theme(panel.background=element_rect(colour=NA,fill="white",size=0.5),
          panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),panel.spacing = unit(0.1,"line"),
          axis.text.x=element_text(color="black"),axis.title.y=element_text(size=12,color="black"), axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color="black"), axis.line=element_line(),axis.ticks.x = element_blank(),
          axis.ticks.length = unit(.15, "cm"),strip.text=element_text(size=10),
          strip.background = element_rect(fill=NA),legend.position = "bottom")
}
  
# BRAF   WT 
# 189  180 

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/figS5c-d.pdf",width = 14,height = 14)
cowplot::plot_grid(plotlist = mutatedGeneBoxplot[-4],ncol = 5,nrow = 4)
dev.off()

braf_boxplot <-  mutatedGeneBoxplot$BRAF
save(braf_boxplot, file= "braf_boxplot.RData")



