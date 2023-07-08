# 
# TCGA metastatic immune-infiltrating 
# 
# REF：Estimating the population abundance of tissue-infiltrating immune and stromal
# cell populations using gene expression
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/03.immuneMechanism/"))
library(ComplexHeatmap)
# MCP score
load( 'd:/Projects/Melanoma_Pyroptosis_new/4.NMF_survival_immune_dataste_基于转移患者显著基因/20220414Res/melabomaPvsM.MCPcounter.rdata')
# PScore
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/PScoreModel_ROC/PyropScore_SKCM.Rdata")
# PScore group
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/heatmap/PScore_group_meta.RData")
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/heatmap/PScore_group_pri.RData")


# ============================================================================ #
# 1. data prepare ####
# ============================================================================ #
data_ES_met <- data_ES_List$metastatic
data_ES_met$sample_id <- rownames(data_ES_met)
table(duplicated(data_ES_met$PRGs_prot))
PScore_group_meta$Group <- paste0("metastatic_",PScore_group_meta$Group)
data_ES_met_group <- merge(data_ES_met[,c("sample_id","PRGs_prot")],PScore_group_meta, by = "PRGs_prot", all.x = T)
data_ES_met_group <- data_ES_met_group[order(data_ES_met_group$PRGs_prot, decreasing = F),]

data_ES_pri <- data_ES_List$primary
data_ES_pri$sample_id <- rownames(data_ES_pri)
table(duplicated(data_ES_pri$PRGs_prot))
PScore_group_pri$Group <- paste0("primary_",PScore_group_pri$Group)
data_ES_pri_group <- merge(data_ES_pri[,c("sample_id","PRGs_prot")],PScore_group_pri, by = "PRGs_prot", all.x = T)
data_ES_pri_group <- data_ES_pri_group[order(data_ES_pri_group$PRGs_prot, decreasing = F),]

heatmap_order <- rbind(data_ES_met_group[,-5],data_ES_pri_group)

# order
cluster_order <- heatmap_order$sample_id
melanomaMCP <- melanomaMCP[,cluster_order]

# scale across cells
melanomaMCP <- t(scale(t(melanomaMCP)))


# ============================================================================ #
# 2.Wilcox test High_PS vs Low_PS ####
# ============================================================================ #
metastatic_High_PS <- PScore_group_meta$barcode[PScore_group_meta$Group=="metastatic_High_PS"]
metastatic_Low_PS <- PScore_group_meta$barcode[PScore_group_meta$Group=="metastatic_Low_PS"]
table(substr(gsub("\\.","_",colnames(melanomaMCP)),1,12) %in% c(metastatic_Low_PS, metastatic_High_PS))

MCP_wilcox <- melanomaMCP[,substr(gsub("\\.","_",colnames(melanomaMCP)),1,12) %in% c(metastatic_Low_PS, metastatic_High_PS)]
group <- ifelse(substr(gsub("\\.","_",colnames(MCP_wilcox)),1,12) %in% metastatic_Low_PS, "Low", "High")
group <- as.factor(group)

wilcoxTestList <- list()
wilcoxTestList[["p.Tcells"]] <- wilcox.test(MCP_wilcox[1,]~group,data = MCP_wilcox)
wilcoxTestList[["p.CD8T"]] <- wilcox.test(MCP_wilcox[2,]~group,data = MCP_wilcox)
wilcoxTestList[["p.CytotoxicLymphocytes"]] <- wilcox.test(MCP_wilcox[3,]~group,data = MCP_wilcox)
wilcoxTestList[["p.Blineage"]] <- wilcox.test(MCP_wilcox[4,]~group,data = MCP_wilcox)
wilcoxTestList[["p.NK"]] <- wilcox.test(MCP_wilcox[5,]~group,data = MCP_wilcox)
wilcoxTestList[["p.Monocytic"]] <- wilcox.test(MCP_wilcox[6,]~group,data = MCP_wilcox)
wilcoxTestList[["p.MyeloidDendriticCells"]] <- wilcox.test(MCP_wilcox[7,]~group,data = MCP_wilcox)
wilcoxTestList[["p.Neutrophils"]] <- wilcox.test(MCP_wilcox[8,]~group,data = MCP_wilcox)
wilcoxTestList[["p.Endothelial"]] <- wilcox.test(MCP_wilcox[9,]~group,data = MCP_wilcox)
wilcoxTestList[["p.Fibroblasts"]] <- wilcox.test(MCP_wilcox[10,]~group,data = MCP_wilcox)

tmp <- list()
for(i in names(wilcoxTestList)){
  if(wilcoxTestList[[i]]$p.value > 0.05){
    tmp[i] <- "ns"
  } else {
    if(wilcoxTestList[[i]]$p.value < 0.001) {
      tmp[i] <- "***"
    } else {
      tmp[i] <- "*"
    }
  }
}
print(tmp)

rownames(melanomaMCP) <- paste0(rownames(melanomaMCP), "\n" ,tmp)



# ============================================================================ #
# 3. heatmap anno prepare ####
# ============================================================================ #

# OS
clinic <- read.delim("d:/Work/Data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
mel_cli <- clinic[clinic$type%in%"SKCM",c("bcr_patient_barcode","OS.time")]  
mel_cli$bcr_patient_barcode <-  gsub("-","_",mel_cli$bcr_patient_barcode)

melanoma <- read.delim("d:/Work/Data/TCGA_preliminary/SKCM/melanoma_mRNA.txt")
mel_cli <- mel_cli[mel_cli$bcr_patient_barcode %in%  gsub("\\.","_",substr(colnames(melanoma),1,12)),]  # 1 person without clinicInfo - Cluster$Sample[Cluster$Sample %in% mel_cli$bcr_patient_barcode ==FALSE]

cluster_order_clean <- substr(gsub("\\.","_",cluster_order),1,12)
cluster_order_clean[cluster_order_clean %in% mel_cli$bcr_patient_barcode ==FALSE]
# "TCGA_GN_A261"

# match
mel_cli <- mel_cli[match(cluster_order_clean,mel_cli$bcr_patient_barcode),]
mel_cli$bcr_patient_barcode[is.na(mel_cli$bcr_patient_barcode)] <- "TCGA_GN_A261"
mel_cli$OS.time <-  as.numeric(mel_cli$OS.time)


# ============================================================================ #
# 4. heatmap ####
# ============================================================================ #
heatmapAnno <- HeatmapAnnotation(
  OS = anno_barplot(mel_cli$OS.time),
  Cluster = heatmap_order$Group,
  col = list(Cluster = c(   # 指定颜色，传入list
    "primary_Low_PS" = "#E18727FF",
    "primary_High_PS" = "#20854EFF", 
    "metastatic_Low_PS" = "#0072B5FF", 
    "metastatic_High_PS" = "#BC3C29FF")),
  annotation_legend_param = list(
    Cluster = list(direction = "horizontal",nrow = 2)))  # 可以通过nrow 和 legend 控制图例方向

p <- Heatmap(melanomaMCP, 
             rect_gp = gpar(col= "white"),
             show_column_names = FALSE,
             cluster_columns = F,
             show_row_dend = F,
             row_names_side = "right",
             column_title = "Microenvironment Cell Populations of Different Clusters in SKCM", 
             name =  "MCPcounter.estimate", 
             top_annotation = heatmapAnno,
             heatmap_legend_param = list(direction = "horizontal"))

pdf("immune_mechanism_TME_tcga_MCP_heatmap.pdf",width = 12,height = 6)
draw(p, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()



