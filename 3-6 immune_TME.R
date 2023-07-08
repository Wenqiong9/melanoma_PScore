# 
# TME
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


# ============================================================================ #
# 1. TME ####
# ============================================================================ #
load("gset_all.add_PScore.Rdata")
TME_subtype1 <- openxlsx::read.xlsx("d:/Work/Data/public_immune/Cancer_Cell_TME4.xlsx",sheet = 1) 

intersect(unique(TME_subtype1$series_id),gset_all$GEO_ID)
# [1] "GSE35640" "GSE19234" "GSE65904" "GSE54467"

# check samples of datasets
for(i in 1:7){
  print(gset_all$GEO_ID[[i]])
  tmp <- gset_all$PScore[[i]]
  tmp$geo_accession <- rownames(tmp)
  if(i == 2){  
    tmp$geo_accession <- substr(gsub("\\.","-",rownames(gset_all$PScore[[i]])),1,12)
  }
  gset_all$PScore[[i]] <- tmp
  print(table(tmp$geo_accession %in% TME_subtype1$Sample))
}
# [1] "inhouse_data"
# FALSE 
# 62 
# [1] "TCGA_metastatic"
# FALSE  TRUE 
# 5   364 
# [1] "GSE19234"
# TRUE 
# 39 
# [1] "GSE35640"
# TRUE 
# 65 
# [1] "GSE54467"
# TRUE 
# 79 
# [1] "GSE65904"
# TRUE 
# 198 
# [1] "PRJEB23709"
# TRUE 
# 72 

# rbind datasets
PScore <- gset_all$PScore  %>%  dplyr::bind_rows() %>% 
  dplyr::select(geo_accession, PScore)  %>% 
  dplyr::mutate(Dataset=rep(gset_all$GEO_ID, unlist(lapply(gset_all$PScore,nrow))))

# merge PScore and TME
TME_PScore <- merge(PScore,TME_subtype1,by.x="geo_accession",by.y="Sample")
TME_PScore <- TME_PScore[,!colnames(TME_PScore)%in%"Patient"]

openxlsx::write.xlsx(TME_PScore,file = "TME_PScore_5_datasets.xlsx",overwrite = T)


# ============================================================================ #
# 2. plot ####
# ============================================================================ #
TME_PScore <- openxlsx::read.xlsx("TME_PScore_5_datasets.xlsx")

table(TME_PScore$MFP)
# D    F   IE IE/F 
# 311  122  204  180 
table(TME_PScore$Dataset)
# GSE19234        GSE35640        GSE54467        GSE65904      PRJEB23709 TCGA_metastatic 
#       39              65              79             198              72             364 


compare <- list(c("D","F"),c("D","IE"),c("D","IE/F"),c("F","IE"),c("F","IE/F"),c("IE","IE/F"))

TME_plotList <- list()
for(i in unique(gset_all$GEO_ID[gset_all$GEO_ID %in% TME_PScore$Dataset])){
  Data <- TME_PScore[TME_PScore$Dataset %in% i,]
  
  TME_plotList[[i]] <- ggplot(data = Data,aes(x=MFP,y=PScore,fill=MFP))+
    geom_boxplot(width = 0.6)+
    scale_x_discrete(limits=c("D","F","IE","IE/F"),labels=paste(c("D","F","IE","IE/F"),"(",table(Data$MFP),")"))+
    scale_fill_manual(limits=c("D","F","IE","IE/F"),values = c("#BC3C29FF","#0072B5FF","#E18727FF","#FFDC91FF"),name="Response",guide="none")+
    stat_compare_means(comparisons = compare,method = "wilcox.test")+
    xlab("") + ylab("Pyroptosis Score")+ labs(title = paste0(i))+
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
}

pdf("immune_mechanism_TME.pdf", width = 15,height = 5)
cowplot::plot_grid(plotlist = TME_plotList,ncol=6)
dev.off()
