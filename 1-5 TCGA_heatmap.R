# 
# aim: heatmap
# annotation: OS, PScore
# patient label: nmf cluster, PScore group
# gene label: NMF feature primary, primary_sig, NMF feature metastatic, metastatic_sig
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
library(ComplexHeatmap)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))

clinic <- read.delim("d:/Projects/Melanoma_Pyroptosis/Genelist/TCGA_ClinicalData_20180420.txt")
load("coxForest/melanoma_expr_list.PRGs_74.RData")  
melanomaExpr <- melanoma_expr_list[["tumor"]]  

# prognostic gene
load("coxForest/uniCox_list.PRG_74.RData")
# NMF clusters
load("nmf_survival/NMF_primary_RPGs74_Res.Rdata")
load("nmf_survival/NMF_metastatic_RPGs74_Res.Rdata")
# PScore
load("PScoreModel_ROC/PyropScore_SKCM.Rdata")
# PScore_group
load("heatmap/PScore_group_meta.RData")  # 1-6
load("heatmap/PScore_group_pri.RData")   # 1-6


# ============================================================================ #
# 1.prepare column annotation ####
# ============================================================================ #

# ============================== #
## 1.1 NMF cluster and PScore ####
# ============================== #
NMF_metastatic <- nmf_metastatic$group_opt  
data_ES_m <- data_ES_zscore_List$metastatic
data_ES_m$Sample <- rownames(data_ES_m)
NMF_metastatic <- merge(NMF_metastatic,data_ES_m[c("PRGs_prot", "Sample")])

NMF_primary <- nmf_primary$group_opt 
data_ES_p <- data_ES_zscore_List$primary
data_ES_p$Sample <- rownames(data_ES_p)
NMF_primary <- merge(NMF_primary,data_ES_p[c("PRGs_prot", "Sample")])

NMF_metastatic$Group <- paste0("M_",NMF_metastatic$Group)
NMF_primary$Group <- paste0("P_",NMF_primary$Group)
Cluster <- rbind(NMF_primary,NMF_metastatic)
Cluster <- Cluster[order(Cluster$Group),]

Cluster$Sample_clean <- gsub("\\.","_",substr(Cluster$Sample,1,12))
Cluster$Sample_clean <- gsub("-","_",substr(Cluster$Sample_clean,1,12))


# ==================== #
## 1.2 PScore group ####
# ==================== #

table(PScore_group_meta$Group)
# Low_PS High_PS 
#    161     192 
table(PScore_group_pri$Group)
# Low_PS High_PS 
#    56      47 

high_samples <- c(PScore_group_meta$barcode[PScore_group_meta$Group == "High_PS"], PScore_group_pri$barcode[PScore_group_pri$Group == "High_PS"])
low_samples <- c(PScore_group_meta$barcode[PScore_group_meta$Group == "Low_PS"], PScore_group_pri$barcode[PScore_group_pri$Group == "Low_PS"])

Cluster$PScore_group[Cluster$Sample_clean %in% high_samples] <- "High_PS"
Cluster$PScore_group[Cluster$Sample_clean %in% low_samples] <- "Low_PS"
table(Cluster$PScore_group)
# High_PS  Low_PS 
#     238     218
intersect(high_samples, low_samples) # duplicated sample
Cluster$PScore_group[Cluster$Sample == "TCGA.ER.A2NF.01A.11R.A18T.07"] <- "High_PS"


# =================================== #
## 1.3 match the order of patients ####
# =================================== #
melanomaExpr <- melanomaExpr[,substr(colnames(melanomaExpr),14,15) %in% c("01","06")]
melanomaExpr <- melanomaExpr[,match(Cluster$Sample,colnames(melanomaExpr))]

mel_cli <- clinic[clinic$type%in%"SKCM",c("bcr_patient_barcode","OS.time")]  
mel_cli$bcr_patient_barcode <-  gsub("-","_",mel_cli$bcr_patient_barcode)
mel_cli <- mel_cli[mel_cli$bcr_patient_barcode %in% gsub("\\.","_",substr(colnames(melanomaExpr),1,12)),] 
mel_cli <- mel_cli[match(Cluster$Sample_clean,mel_cli$bcr_patient_barcode),]
mel_cli$OS.time <-  as.numeric(mel_cli$OS.time)


# ============================ #
## 1.4 scale across samples ####
# ============================ #
melanomaExpr_metastatic <- melanomaExpr[,1:361]
melanomaExpr_metastatic <- t(scale(t(melanomaExpr_metastatic)))

melanomaExpr_primary <- melanomaExpr[,362:464]
melanomaExpr_primary <- t(scale(t(melanomaExpr_primary)))


# ============================================================================ #
# 2.heatmap-metastatic ####
# ============================================================================ #

# row anno
NMF_feature <- ifelse(rownames(melanomaExpr) %in% nmf_metastatic$feat_opt, "Selected", "Not selected")
NMF_feature <- factor(NMF_feature,levels = c("Selected", "Not selected"))
Univariate_COX <- ifelse(rownames(melanomaExpr_metastatic) %in% uniCox_list$metastatic$gene[uniCox_list$metastatic$Risk_group=="Protect"], "Protect", 
                         ifelse(rownames(melanomaExpr_metastatic) %in% uniCox_list$metastatic$gene[uniCox_list$metastatic$Risk_group=="Risk"], "Risk", "Not significant"))
Univariate_COX <- factor(Univariate_COX,levels = c("Risk", "Protect", "Not significant"))

heatmapAnnoRow <- rowAnnotation(
  NMF_feature = NMF_feature,
  Univariate_COX = Univariate_COX,
  col = list(NMF_feature = c("Selected" = "#FFDC91FF", "Not selected" = "#FFDC9166"),
             Univariate_COX  = c("Protect" = "#0072B5EE", "Risk" = "#BC3C29EE", "Not significant" = "#0072B599")))


heatmapAnno <- HeatmapAnnotation(
  OS = anno_barplot(mel_cli[1:361,]$OS.time),
  PScore = Cluster[1:361,]$PRGs_prot,
  PScore_group = Cluster[1:361,]$PScore_group,
  NMF_cluster = Cluster[1:361,]$Group,
  col = list(NMF_cluster = c("M_1" = "#BC3C29FF", "M_2" = "#0072B5FF", "M_3" = "#E18727FF"),
             PScore_group = c("High_PS"="#20854EFF", "Low_PS"="#20854E99"),
             PScore = circlize::colorRamp2(c(-2,0,3), c("blue", "white", "red"))))

Heatmap(melanomaExpr_metastatic, 
        rect_gp = gpar(col= "#FFFFFF00"),
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_dend = F,
        row_names_side = "left",
        row_km = 3,
        column_title="Expression of Pyroptpsis-related Genes",
        name =  "Expressions",
        top_annotation = heatmapAnno,
        left_annotation = heatmapAnnoRow)


# ============================================================================ #
# 3.heatmap-primary ####
# ============================================================================ #

# row anno
NMF_feature <- ifelse(rownames(melanomaExpr) %in% nmf_primary$feat_opt, "Selected", "Not selected")
NMF_feature <- factor(NMF_feature,levels = c("Selected", "Not selected"))
Univariate_COX <- ifelse(rownames(melanomaExpr_primary) %in% uniCox_list$primary$gene[uniCox_list$primary$Risk_group=="Protect"], "Protect", "Not significant")
Univariate_COX <- factor(Univariate_COX,levels = c("Protect", "Not significant"))

heatmapAnnoRow <- rowAnnotation(
  NMF_feature = NMF_feature,
  Univariate_COX = Univariate_COX,
  col = list(NMF_feature = c("Selected" = "#FFDC91FF", "Not selected" = "#FFDC9166"),  
             Univariate_COX  = c("Protect" = "#0072B5EE", "Not significant" = "#0072B599")))

heatmapAnno <- HeatmapAnnotation(
  OS = anno_barplot(mel_cli[362:464,]$OS.time),
  PScore = Cluster[362:464,]$PRGs_prot,
  PScore_group = Cluster[362:464,]$PScore_group,
  NMF_cluster = Cluster[362:464,]$Group,
  col = list(NMF_cluster = c("P_1" = "#BC3C29FF", "P_2" = "#0072B5FF", "P_3" = "#E18727FF"),
             PScore_group = c("High_PS"="#20854EFF", "Low_PS"="#20854E99"),
             PScore = circlize::colorRamp2(c(-2,0,3), c("blue", "white", "red"))))

Heatmap(melanomaExpr_primary, 
        rect_gp = gpar(col= "#FFFFFF00"), 
        cluster_columns = FALSE,
        show_column_names = FALSE, 
        show_row_dend = F, 
        row_names_side = "left",
        row_km = 3,
        column_title="Expression of Pyroptpsis-related Genes",
        name =  "Expressions",
        top_annotation = heatmapAnno,  
        left_annotation = heatmapAnnoRow)
