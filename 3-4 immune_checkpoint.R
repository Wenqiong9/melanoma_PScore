# 
# Figure 4 immune checkpoints
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
# 1. data prepare ####
# ============================================================================ #
load("gset_all.add_PScore.Rdata")

# InAC <- read.csv("d:/Data/public_immune/immuneGene_Signatures/InAcMarker_New.csv")
# InAC <- read.csv("/work/yye/AliRstudio/Data/Public/Immune/ImmuneGeneSigantures/InAcMarker_New.csv")

InAC <- openxlsx::read.xlsx("d:/Work/Data/public_immune/immuneGene_Signatures/Immunecheckpoint_The Immune Landscape of Cancer (TCGA).xlsx")
InAC$Immune.Checkpoint[InAC$Immune.Checkpoint %in% c("Stumulatory","Stimulaotry")] <- "Stimulatory"
InAC <- InAC[InAC$Immune.Checkpoint %in% "Stimulatory",]
colnames(InAC)[colnames(InAC)=="HGNC.Symbol"] <- "marker"


for(i in 1:7){
  expr_tmp <- gset_all$gset_data_matrix_symbol[[i]]
  print(gset_all$GEO_ID[[i]])
  if("CD27" %in% rownames(expr_tmp) == FALSE){rownames(expr_tmp)[rownames(expr_tmp) == "TNFRSF7"] <- "CD27"}
  if("CD70" %in% rownames(expr_tmp) == FALSE){rownames(expr_tmp)[rownames(expr_tmp) == "TNFSF7"] <- "CD70"}
  print(InAC$marker[InAC$marker %in% rownames(expr_tmp)==FALSE])
}

# expr_tmp <- gset_all$gset_data_matrix_symbol[[5]]
# c("CD27", "S152", "S152. LPFS2", "T14", "TNFRSF7","Tp55") %in% rownames(expr_tmp)
# c("CD70", "CD27-L", "CD27L", "CD27LG", "LPFS3", "TNFSF7", "TNLG8A")%in% rownames(expr_tmp)
# c("IFNA1", "IFL", "IFN", "IFN-ALPHA", "IFN-alphaD", "IFNA13", "IFNA@", "leIF D")%in% rownames(expr_tmp)
# c("IFNA2", "IFN-alpha-2", "IFN-alphaA", "IFNA", "IFNA2B", "leIF A")%in% rownames(expr_tmp)
# c("IL2", "IL-2", "TCGF", "lymphokine")%in% rownames(expr_tmp)

# [1] "inhouse_data"
# character(0)
# [1] "TCGA_metastatic"
# character(0)
# [1] "GSE19234"
# character(0)
# [1] "GSE35640"
# character(0)
# [1] "GSE54467"
# [1] "IFNA1" "IFNA2" "IL2"  
# [1] "GSE65904"
# character(0)
# [1] "PRJEB23709"
# character(0)

for(i in 1:7){
  print(gset_all$GEO_ID[[i]])
  tmp <- gset_all$PScore[[i]]
  tmp$geo_accession <- rownames(tmp)
  gset_all$PScore[[i]] <- tmp
}

lapply(gset_all$PScore, head)


# ============================================================================ #
# 2. analysis ####
# ============================================================================ #

# merge and cor

InAc_PS_Corr <- purrr::map2(.x=gset_all$gset_data_matrix_symbol,.y=gset_all$PScore,function(.x,.y){
  
  Data  <- .x %>% data.frame()
  Data2 <- .y %>% data.frame()
  Data  <- Data[rownames(Data)%in%c(InAC$marker),] %>% t() %>% data.frame() %>% tibble::rownames_to_column(var="geo_accession")
  Data_m <- merge(Data,Data2,by="geo_accession")
  
  Cor <- data.frame(Gene=colnames(Data)[-1])
  Cor[,c("estimate.rho","pvalue")] <- t(apply(Data_m[,Cor$Gene], 2, function(x)unlist(my.cor.test(as.numeric(x),as.numeric(Data_m$PScore),method="pearson"))))   # spearman
  Cor$FDR <- p.adjust(Cor$pvalue,method = "fdr")
  return(Cor)
  
})
# Warning message:
#   In cor(x, y) : 标准差为零
# tmp=gset_all$gset_data_matrix_symbol[[7]];.y=gset_all$PScore[[7]]
# tmp <- tmp[rownames(tmp)%in%c(InAC$marker),]
# the reason of warning is the expr of TNF is 0.

InAc_PS_Corr_5_Datasets <- InAc_PS_Corr %>% dplyr::bind_rows() %>% dplyr::mutate(Dataset=rep(gset_all$GEO_ID,times=unlist(lapply(InAc_PS_Corr,nrow))))
write.csv(InAc_PS_Corr_5_Datasets,file = "ImmuneCheckPoint_PS_Corr_5_Datasets.csv",quote = F)


# ============================================================================ #
# 3. plot ####
# ============================================================================ #
InAc_PS_Corr_5_Datasets <- read.csv("ImmuneCheckPoint_PS_Corr_5_Datasets.csv")

InAc_PS_Corr_5_Datasets$FDR_log10 <- -log10(InAc_PS_Corr_5_Datasets$FDR + 10^-10)
InAc_PS_Corr_5_Datasets$FDR_log10[InAc_PS_Corr_5_Datasets$FDR_log10 > -log(0.05)] <- -log10(1e-10)

InAc_PS_Corr_5_Datasets$class <- InAC$class[match(InAc_PS_Corr_5_Datasets$Gene,InAC$marker)]
y_order <- t(sapply(split(InAc_PS_Corr_5_Datasets[,c("estimate.rho","FDR_log10")],InAc_PS_Corr_5_Datasets$Gene),function(x) 
  c(nrow(x[x$estimate.rho > 0 & x$FDR_log10 > -log10(0.05),]),nrow(x[x$estimate.rho < 0 & x$FDR_log10 > -log10(0.05),])))) %>% 
  data.frame() %>% 
  set_colnames(c("Pos","Neg"))
y_order <- y_order[order(y_order$Pos-y_order$Neg),]

estimate_range <- max(abs(range(InAc_PS_Corr_5_Datasets$estimate.rho,na.rm = T)))

pdf("immune_mechanism_checkpoint.pdf",width=9,height = 5) 
ggplot(InAc_PS_Corr_5_Datasets,aes(y=Dataset,x=Gene,color=estimate.rho,fill=estimate.rho))+
  geom_point(aes(size=FDR_log10),shape=23)+
  scale_color_gradientn(limit=c(-estimate_range,estimate_range),colors= colorRampPalette(c("#0072B5FF","white","red"),space="rgb")(100),breaks=seq(-1,1,length.out = 5),labels=c("-1","-0.5","0","0.5","1"),name="Rs")+
  scale_fill_gradientn(limit=c(-estimate_range,estimate_range),colors= colorRampPalette(c("#0072B5FF","white","red"),space="rgb")(100),breaks=seq(-1,1,length.out = 5),labels=c("-1","-0.5","0","0.5","1"),name="Rs")+
  scale_size_continuous(limit=c(0,10),range = c(0.1, 5),breaks=c(-log10(0.05),5,10),labels=c("0.05","1e-5","<1e-10"),name="FDR")+
  geom_point(data=InAc_PS_Corr_5_Datasets[InAc_PS_Corr_5_Datasets$FDR_log10 > 1.3,],aes(y=Dataset,x=Gene,size=FDR_log10),shape=23,color="black") +
  coord_fixed()+
  scale_y_discrete(limit=rev(names(table(InAc_PS_Corr_5_Datasets$Dataset))))+
  scale_x_discrete(limit=rownames(y_order))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(size=8,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),legend.position = "bottom",
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()

