# Aim: Correlation of expression of PRGs 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))

library(corrplot)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# ============================================================================ #
# 1.correlation of 74 gene expression####
# ============================================================================ #
load("coxForest/Expr_Survival_list.PRG_74.RData")
names(Expr_Survival_list)<- c("SKCM (Primary)","SKCM (Metastastic)")

p_value <- list()
exprCor <- list()
for(j in names(Expr_Survival_list)){
  mydata <- Expr_Survival_list[[j]]
  mydata <- mydata[,-c(1,2)]; head(mydata); print(dim(mydata))
  
  mydata_cor <- cor(mydata)
  p.mat <- cor.mtest(mydata_cor); head(p.mat[, 1:5])
  p_value[[j]] <- p.mat
  exprCor[[j]] <- mydata_cor
  
  # pdf(paste0("geneExprCorrelation/PRGs_74.", j, ".pdf"), width = 7, height = 8)
  corrplot(mydata_cor, method = "color", mar = c(0, 0, 2, 0),
           title = paste0(j,"\nHclust.method: ward.D"),
           tl.col = "black", tl.cex = 0.6,  number.cex = 0.3, tl.srt = 45,
           order = "hclust", addrect = 3, hclust.method = "ward.D", 
           p.mat=p.mat,insig = "label_sig",
           sig.level = c(0.001, 0.01, 0.05), pch.cex = .3, pch.col = "black")
  # dev.off()
}
# [1] 103  74
# [1] 361  74

openxlsx::write.xlsx(p_value, row.names = TRUE, file = "geneExprCorrelation/PRGs_74.p_value.xlsx")
openxlsx::write.xlsx(exprCor, row.names = TRUE, file = "geneExprCorrelation/PRGs_74.exprCor.xlsx")


# ============================================================================ #
# 2.correlation of 74 gene expression####
# ============================================================================ #
load("coxForest/uniCoxSig_list.RData")
meta_sig_genes <- uniCoxSig_list$metastatic$gene


p_value <- list()
exprCor <- list()
corrplot_list <- list()
for(j in names(Expr_Survival_list)){
  mydata <- Expr_Survival_list[[j]]
  mydata <- mydata[,-c(1,2)]
  mydata <- mydata[,meta_sig_genes]; head(mydata); print(dim(mydata))
  
  mydata_cor <- cor(mydata)
  p.mat <- cor.mtest(mydata_cor); head(p.mat[, 1:3])
  # p_value[[j]] <- p.mat
  # exprCor[[j]] <- mydata_cor
  
  pdf(paste0("geneExprCorrelation/PRGs_meta_Sig34.", j, ".pdf"), width = 6, height = 5)
  corrplot_list[[j]] <- corrplot(mydata_cor, method = "color", mar=c(0, 0, 2, 0),
           title = paste0(j,"\nHclust.method: ward.D"),
           tl.col = "black", tl.cex = 0.6,  number.cex = 0.3, tl.srt = 45,
           order = "hclust", addrect = 4, hclust.method = "ward.D", 
           p.mat = p.mat, insig = "label_sig",
           sig.level = c(0.001, 0.01, 0.05), pch.cex = .6, pch.col = "black")
  dev.off()

}
