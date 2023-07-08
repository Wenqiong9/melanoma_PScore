# heatmap 
# 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
library(ggplot2)
library(ggpubr)
library(maftools)
library(ComplexHeatmap)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/02.mutation/"))

load("D:/Work/Data/TCGAmutations/extdata/MC3/SKCM.RData")
SKCM <- tcga_skcm_mc3@data

load(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/PScoreModel_ROC/PyropScore_SKCM.Rdata"))

load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")
genelist <- list(PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])
candidate_genes <- genelist$PRGs_prot


# ============================================================================ #
# oncoplot: top10 and PRGs ####
# ============================================================================ #
skcm_metastatic_Maf <- subsetMaf(maf = tcga_skcm_mc3, tsb = substr(gsub("\\.", "-", rownames(data_ES_List$metastatic)),1,12)) 
candidate_genes[candidate_genes %in% skcm_metastatic_Maf@data$Hugo_Symbol==FALSE]

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/figS5a-b_oncoplot.pdf",width = 7, height = 7)
maftools::oncoplot(maf = skcm_metastatic_Maf, top = 10)
maftools::oncoplot(maf = skcm_metastatic_Maf, genes = candidate_genes)
dev.off()


# ============================================================================ #
# heatmap ####
# ============================================================================ #
top10Gene <- names(table(skcm_metastatic_Maf@data$Hugo_Symbol)[order(table(skcm_metastatic_Maf@data$Hugo_Symbol), decreasing = T)][c(1:6,8:10,29)])
# [1] "TTN"     "MUC16"   "DNAH5"   "PCLO"    "LRP1B"   "GPR98"   "DNAH7"   "PKHD1L1" "CSMD1"   "BRAF"  


load("mutation_complexHeatmap.RData") # data_ES_m
colnames(data_ES_m)[colnames(data_ES_m)=="BRAF_status"] <- "BRAF"

mut_total <- substr(unique(skcm_metastatic_Maf@data$Tumor_Sample_Barcode_full),1,16)
for(i in c(top10Gene, "PTEN","NRAS","KIT","TERT",'NLRP1', 'NLRP3', "NLRP7","TLR4")){
  # i = "TTN"
  SKCM_mut <- skcm_metastatic_Maf@data[skcm_metastatic_Maf@data$Hugo_Symbol == i,]
  mut_candidate <- substr(unique(SKCM_mut$Tumor_Sample_Barcode_full),1,16)
  data_ES_m[[i]] <- ifelse(data_ES_m$barcode %in% mut_candidate, i, ifelse(data_ES_m$barcode %in% mut_total, "WT", NA))
  tmp <- table(data_ES_m[[i]])
  print(paste(i, tmp[1]/364))
}

data_ES_m <- data_ES_m[order(data_ES_m$PRGs_prot),]
head(data_ES_m)


heatmapAnno <- HeatmapAnnotation(
  TTN = data_ES_m$TTN, MUC16 = data_ES_m$MUC16, DNAH5 = data_ES_m$DNAH5,
  PCLO = data_ES_m$PCLO, LRP1B = data_ES_m$LRP1B, GPR98 = data_ES_m$GPR98,
  PKHD1L1 = data_ES_m$PKHD1L1, CSMD1 = data_ES_m$CSMD1, DNAH7 = data_ES_m$DNAH7,
  NLRP1 = data_ES_m$NLRP1, NLRP3 = data_ES_m$NLRP3, NLRP7 = data_ES_m$NLRP7, TLR4 = data_ES_m$TLR4,
  PTEN = data_ES_m$PTEN,TERT = data_ES_m$TERT,KIT = data_ES_m$KIT,NRAS = data_ES_m$NRAS,
  BRAF = data_ES_m$BRAF, PScore_group = data_ES_m$group,
  col = list(PScore_group = c("High_PS"="#BC3C29FF", "Low_PS"="#0072B5FF"),
             TTN = c("TTN" = "#7876BEFF", "WT" = "#7876BE33"), MUC16 = c("MUC16" = "#6F99ADFF", "WT" = "#6F99AD33"), 
             DNAH5 = c("DNAH5" = "#FFDC91FF", "WT" = "#FFDC9133"), PCLO = c("PCLO" = "#3B3B3BFF", "WT" = "#3B3B3B33"), 
             LRP1B = c("LRP1B" = "#E64B35FF", "WT" = "#E64B3533"), GPR98 = c("GPR98" = "#EE4C97FF", "WT" = "#EE4C9733"),
             PKHD1L1 = c("PKHD1L1" = "#003C67FF", "WT" = "#003C6733"), CSMD1 = c("CSMD1" = "#7AA6DCFF", "WT" = "#7AA6DC33"), 
             DNAH7 = c("DNAH7" = "#00A087FF", "WT" = "#00A08733"), NLRP1 = c("NLRP1" = "#868686FF", "WT" = "#86868633"),
             NLRP3 = c("NLRP3" = "#8F7700FF", "WT" = "#8F770033"), NLRP7 = c("NLRP7" = "#EFC000FF", "WT" = "#EFC00033"),
             TLR4 = c("TLR4" = "darkgrey", "WT" = "lightgrey"),
             PTEN =  c("PTEN" = "#BC3C29FF", "WT" = "#BC3C2933"), NRAS =  c("NRAS" = "#20854EFF", "WT" = "#20854E33"),
             KIT =  c("KIT" = "#E18727FF", "WT" = "#E1872733"), TERT =  c("TERT" = "#0072B5FF", "WT" = "#0072B533"),
             BRAF = c("BRAF" = "#A73030FF", "WT" = "#A7303033")),
  annotation_name_side = "left"
)

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig2a_heatmap.pdf",width = 16,height = 4.6)
Heatmap(t(as.matrix(data_ES_m$PRGs_prot)), 
        rect_gp = gpar(col= "#FFFFFF00"), 
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_dend = F, 
        row_names_side = "left",
        column_title="Mutation of Metastatic Melanoma",
        name = "Pyroptosis_Score", 
        top_annotation = heatmapAnno,
        left_annotation = NULL)
dev.off()

