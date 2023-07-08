# 
# FIGURE 4 GEP CYT
# 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/03.immuneMechanism/"))
library(tibble)
library(ggplot2)
library(ggpubr)
library(magrittr)
my.cor.test <- function(...) {
  obj<-try(cor.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}


# ============================================================================ #
# 1.data prepare ####
# ============================================================================ #
load("gset_all.add_PScore.Rdata")

GEP_list <- read.delim("d:/Work/Data/public_immune/immuneGene_Signatures/GEP_List.txt")
GEP_list <- list(GEP=GEP_list$Gene)
CYT <- list(CYT=c("GZMA","PRF1"))

for(i in 1:7){
  expr_tmp <- gset_all$gset_data_matrix_symbol[[i]]
  if("CD27" %in% rownames(expr_tmp) == FALSE){rownames(expr_tmp)[rownames(expr_tmp) == "TNFRSF7"] <- "CD27"}
  if("IDO1" %in% rownames(expr_tmp) == FALSE){rownames(expr_tmp)[rownames(expr_tmp) == "INDO"] <- "IDO1"}
  if("TIGIT" %in% rownames(expr_tmp) == FALSE){rownames(expr_tmp)[rownames(expr_tmp) == "VSIG9"] <- "TIGIT"}
  print(paste0(gset_all$GEO_ID[[i]],": ",c(GEP_list$GEP,"GZMA","PRF1")[c(GEP_list$GEP,"GZMA","PRF1") %in% rownames(expr_tmp)==FALSE]))
}
# expr_tmp <- gset_all$gset_data_matrix_symbol[[3]]
# c("HLA-DRB1", "DRB1", "HLA-DR1B", "HLA-DRB", "SS1") %in% rownames(expr_tmp)

# [1] "inhouse_data: "
# [1] "TCGA_metastatic: "
# [1] "GSE35640: HLA-DRB1"
# [1] "GSE53118: "
# [1] "GSE54467: "
# [1] "GSE65904: "
# [1] "PRJEB23709: "

for(i in 1:7){
  print(gset_all$GEO_ID[[i]])
  tmp <- gset_all$PScore[[i]]
  tmp$geo_accession <- rownames(tmp)
  gset_all$PScore[[i]] <- tmp
}


# ============================================================================ #
# 2.GEP CYT analysis ####
# ============================================================================ #
# GEP
gset_all <- gset_all %>% dplyr::mutate(GEP=purrr::map(.x=gset_data_matrix_symbol,function(.x){
  Data <- .x %>% as.matrix()
  data_ES <- GSVA::gsva(Data,GEP_list)
  # data_ES <- data_ES %>% t() %>% data.frame()
  # data_ES$geo_accession <- rownames(data_ES)
  return(data_ES)
}))

# CYT
gset_all <- gset_all %>% dplyr::mutate(CYT=purrr::map(.x=gset_data_matrix_symbol,function(.x){
  Data <- .x %>% data.frame()   # note: - -> . colnames 
  Data_F <- Data[c("GZMA","PRF1"),]
  tempy <- apply(Data_F, 2, function(x){
    rt=sqrt(x[1]*x[2])
  })
  tempy_ES <- data.frame(t(data.frame(tempy)))
  rownames(tempy_ES) <- "CYT"
  tempy_ES <- tempy_ES %>% t() %>% data.frame()
  tempy_ES$geo_accession <- rownames(tempy_ES)
  return(tempy_ES)
}))

# GEP_CYT
gset_all <- gset_all %>% dplyr::mutate(GEP_CYT=purrr::map2(.x=GEP,.y=CYT,function(.x,.y){
  GEP <- .x %>% t() %>% data.frame()
  GEP$geo_accession <- rownames(GEP)
  CYT <- .y %>% data.frame()
  GEP_CYT <- merge(GEP,CYT,by="geo_accession")
  return(GEP_CYT)
}))

# GEP_CYT and PS
gset_all <- gset_all %>% dplyr::mutate(GEP_CYT=purrr::map2(.x=GEP_CYT,.y=PScore  ,function(.x,.y){
  GEP_CYT <- .x %>% data.frame()
  .y <- .y %>% data.frame()
  GEP_CYT <- merge(GEP_CYT,.y[c("PScore","geo_accession")],by = "geo_accession")
  return(GEP_CYT)
}))

# corr of GEP_CYT and PS
GEP_CYT_PS_Corr <- purrr::map(.x=gset_all$GEP_CYT,function(.x){
  Data <- .x %>% data.frame()
  Cor <- data.frame(Gene=c("CYT","GEP"))
  Cor[,c("estimate.rho","pvalue")] <- t(apply(Data[c("CYT","GEP")], 2, function(x)unlist(my.cor.test(as.numeric(x),as.numeric(Data$PScore),method="pearson"))))
  return(Cor)
})

GEP_CYT_PS_Corr_5_Datasets <- GEP_CYT_PS_Corr %>% dplyr::bind_rows() %>% dplyr::mutate(Dataset=rep(gset_all$GEO_ID,each=2))
write.csv(GEP_CYT_PS_Corr_5_Datasets,file="GEP_CYT_PS_Corr_5_Datasets.csv",quote = F)  # 7 datasets in fact


# ============================================================================ #
# 3. plot ####
# ============================================================================ #
GEP_CYT_PS_Corr_5_Datasets <- read.csv("GEP_CYT_PS_Corr_5_Datasets.csv")

GEP_CYT_PS_Corr_5_Datasets$log10P <- -log10(GEP_CYT_PS_Corr_5_Datasets$pvalue)

x_order2 <- lapply(split(GEP_CYT_PS_Corr_5_Datasets[,c("estimate.rho")],GEP_CYT_PS_Corr_5_Datasets$Dataset), mean) %>% 
  unlist() %>% sort() %>% names()

pdf("immune_mechanism_GEP_CYT.pdf",width = 8,height= 6)
ggplot(GEP_CYT_PS_Corr_5_Datasets,aes(x=estimate.rho,y = Dataset,fill=log10P))+ 
  geom_bar(width=0.5,stat = "identity")+
  scale_y_discrete(limits = x_order2)+
  facet_wrap(~Gene,scales = "free_y",nrow = 1)+
  scale_fill_gradient2(low = "blue",mid="purple", high = "red",midpoint = -log10(0.05),name="-log10 (P)")+
  ylab("")+xlab("")+
  theme(panel.background=element_rect(colour=NA,fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.x=element_text(size=10,color="black",vjust=0.5,hjust=0.5),
        axis.title.y=element_text(size=10,color="black",vjust=0.5,hjust=0.5),
        axis.text.x = element_text(angle = 90,size=10,color="black",vjust=0.5,hjust=1),
        #axis.text.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        legend.position = "right")
dev.off()
