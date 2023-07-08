# single cell


rm(list = ls())
setwd("d:/Work/Projects/BMC_Revised/Results/04.singleCell_ST/")

library(Seurat)
library(ggplot2)
library(dplyr)
library(MySeuratWrappers) 

SKCM_TISHdata <- readRDS('d:/Work/Data/public_scRNA/TISH_DB_Melanoma/TISH_DB_ScRNAseqData.rds.gz')
SKCM_TISHdata <- SKCM_TISHdata[c(1,6),]
singleCell_PScore <- readRDS("d:/Work/Projects/BMC_Revised/Results/04.singleCell_ST/singleCell_PScore_GSE115978_GSE75260.rds")

load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")
markers <- uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"]
markers[markers == "DFNB59"] <- "PJVK"  


# ==============================================================================
# 1.data prepare ####
# ==============================================================================

SKCM_TISHdata_PScore <- SKCM_TISHdata

for(i in 1:2){
  GSE_PScore <- singleCell_PScore$Pyroptosis_Score[[i]]
  GSE_seurat <- SKCM_TISHdata$SeuratObject[[i]]
  
  GSE_meta <- GSE_seurat@meta.data
  print(table(rownames(GSE_meta) == rownames(GSE_PScore)))
  GSE_meta <- cbind(GSE_meta,GSE_PScore)   
  GSE_seurat@meta.data <- GSE_meta
  
  PScoreMatrix <- apply(GSE_PScore,2,function(x)signif(x,digits = 3))
  PScoreMatrix <- t(as.matrix(PScoreMatrix))
  PScoreMatrix_Seurat <- CreateSeuratObject(counts = PScoreMatrix )

  PScoreMatrix_Seurat@meta.data <- GSE_seurat@meta.data   
  GSE_seurat@assays$Pyroptosis_Score <- PScoreMatrix_Seurat@assays$RNA  
  
  SKCM_TISHdata_PScore$SeuratObject[[i]] <- GSE_seurat
}

rm(SKCM_TISHdata,singleCell_PScore,GSE_meta,GSE_PScore,GSE_seurat,PScoreMatrix,PScoreMatrix_Seurat)

GSE72056_plotlist <- list()
GSE115978_plotlist <- list()


# ==============================================================================
# 2.GSE115978 ####
# ==============================================================================

GSE_seurat <- SKCM_TISHdata_PScore$SeuratObject[[1]] 
GSE_id <- gsub("/.*$","",SKCM_TISHdata_PScore$MetaDataLists[[1]])

Idents(GSE_seurat) <- GSE_seurat$Stage
GSE_seurat <- subset(GSE_seurat,idents=c("Metastatic")) 
Idents(GSE_seurat) <- GSE_seurat$Treatment
GSE_seurat <- subset(GSE_seurat,idents=c("None")) 

Idents(GSE_seurat) <- GSE_seurat$Celltype..minor.lineage.
table(GSE_seurat$Celltype..minor.lineage.)
#   B       CD4Tn      CD8Tex Endothelial Fibroblasts   Malignant    Monocyte          NK         pDC     Tprolif 
# 424         381         698          36          38        1107          75          45          25          76 


GSE115978_plotlist[["DimPlot"]] <-
DimPlot(GSE_seurat,reduction = 'umap',cols = paletteer::paletteer_d("ggsci::nrc_npg")[c(4,5,7,3,1,2,10,9,6,8)],
        pt.size = 1,label = T,repel = T)+
  ggtitle(GSE_id)+   
  theme(plot.title = element_text(hjust = 0,size=12), 
        axis.title= element_text(size=10))+ 
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line = element_blank())


GSE_meta <- GSE_seurat@meta.data
x_order <- lapply(split(GSE_meta$PScore,GSE_meta$Celltype..minor.lineage.), median) %>% unlist()
x_order <- names(x_order)[order(x_order)]

GSE115978_plotlist[["VlnPlot"]] <- 
VlnPlot(GSE_seurat,features =c("PScore"),pt.size = 0.001,
        cols = paletteer::paletteer_d("ggsci::nrc_npg")[c(4,5,7,3,1,2,10,9,6,8)]) +
  labs(title=paste0(GSE_id,"_Pyropotosis_Score"))+
  NoLegend()+
  scale_x_discrete(limits=x_order)


GSE115978_plotlist[["VlnPlot_markers"]] <-
VlnPlot(GSE_seurat, features = markers, sort = TRUE,
        stacked=T,pt.size= 0, x.lab = '', y.lab = '', 
        cols =paletteer::paletteer_d("ggsci::nrc_npg")[c(4,5,7,3,1,2,10,9,6,8)],
        direction = "horizontal")+ 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

  
GSE115978_plotlist[["FeaturePlot"]] <-
FeaturePlot(GSE_seurat,reduction = "umap",features =c("PScore"),
            cols = c("white","lightgrey", 'blue', "seagreen2"),
            pt.size = 1,label=F)+ 
  labs(title=paste0(GSE_id,"_PScore"))+
  theme(plot.title = element_text(hjust = 0,size=12),
        axis.title = element_text(size=10))+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank())


GSE115978_plotlist[["DotPlot"]] <-
DotPlot(GSE_seurat,features = c("PScore"),cols =c("RdBu"))+ 
  NoLegend()+
  theme(axis.text.x  = element_text(angle = 0,hjust = 0.5),
        legend.position = 'right',legend.text = element_text(size = 6),
        axis.title = element_blank(),
        title = element_text(size = 12))+
    scale_y_discrete(limits=rev(x_order))+
    scale_x_discrete(label="Pyroptosis_Score")+
    labs(title=GSE_id)



# ==============================================================================
# 3. GSE72056 ####
# ==============================================================================

GSE_seurat <- SKCM_TISHdata_PScore$SeuratObject[[2]] 
GSE_id <- gsub("/.*$","",SKCM_TISHdata_PScore$MetaDataLists[[2]])

Idents(GSE_seurat) <- GSE_seurat$Celltype..minor.lineage.
table(GSE_seurat$Celltype..minor.lineage.)
#   B       CD4Tn      CD8Tex Endothelial Fibroblasts   Malignant    Monocyte         Tfh         Th1     Tprolif 
# 642         604        1002          91          87        1365         221         457          84          92 


GSE72056_plotlist[["DimPlot"]] <-
  DimPlot(GSE_seurat,reduction = 'umap',cols = c("#3C5488FF", "#692294", "#F39B7FFF", "#FFE8FF",paletteer::paletteer_d("ggsci::nrc_npg")[c(7,3,1,2,10,8,6,9)]),
          pt.size = 1,label = T,repel = T)+
  ggtitle(GSE_id)+
  theme(plot.title = element_text(hjust = 0,size=12),
        axis.title= element_text(size=10))+ 
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line = element_blank()) 


GSE_meta <- GSE_seurat@meta.data
x_order <- lapply(split(GSE_meta$PScore,GSE_meta$Celltype..minor.lineage.), median) %>% unlist()
x_order <- names(x_order)[order(x_order)]

GSE72056_plotlist[["VlnPlot"]] <- 
  VlnPlot(GSE_seurat,features =c("PScore"),pt.size = 0.001,
          cols =  c("#3C5488FF", "#692294", "#F39B7FFF", "#FFE8FF", paletteer::paletteer_d("ggsci::nrc_npg")[c(7,3,1,2,10,8,6,9)])) +
  labs(title=paste0(GSE_id,"_Pyropotosis_Score"))+
  NoLegend()+
  scale_x_discrete(limits=x_order)


GSE72056_plotlist[["VlnPlot_markers"]] <-
  VlnPlot(GSE_seurat, features = markers, sort = TRUE,
          stacked=T,pt.size= 0, x.lab = '', y.lab = '',
          cols = c("#3C5488FF", "#692294", "#F39B7FFF", "#FFE8FF", paletteer::paletteer_d("ggsci::nrc_npg")[c(7,3,1,2,10,9)]),       # 颜色 
          direction = "horizontal")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())


GSE72056_plotlist[["FeaturePlot"]] <-
  FeaturePlot(GSE_seurat,reduction = "umap",features =c("PScore"),
              cols = c("white","lightgrey", 'blue', "seagreen2"),
              pt.size = 1,label=F)+
  labs(title=paste0(GSE_id,"_PScore"))+
  theme(plot.title = element_text(hjust = 0,size=12),
        axis.title = element_text(size=10))+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank())


GSE72056_plotlist[["DotPlot"]] <-
  DotPlot(GSE_seurat,features = c("PScore"),cols =c("RdBu"))+ #,split.by = "Treatment"
  theme(axis.text.x  = element_text(angle = 0,hjust = 0.5),
        legend.position = 'right',legend.text = element_text(size = 6),
        axis.title = element_blank(),
        title = element_text(size = 12))+
  scale_y_discrete(limits=rev(x_order))+
  scale_x_discrete(label="Pyroptosis_Score")+
  labs(title=GSE_id)



# ==============================================================================
# 4.pdf ####
# ==============================================================================
p0 <- cowplot::plot_grid(plotlist = list(GSE115978_plotlist$DimPlot,GSE72056_plotlist$DimPlot), ncol=2)
p1 <- cowplot::plot_grid(plotlist = list(GSE115978_plotlist$FeaturePlot,GSE72056_plotlist$FeaturePlot), ncol=2)
p2 <- cowplot::plot_grid(plotlist = list(GSE115978_plotlist$VlnPlot,GSE72056_plotlist$VlnPlot), ncol=2)
p3 <- cowplot::plot_grid(plotlist = list(GSE115978_plotlist$DotPlot,GSE72056_plotlist$DotPlot), ncol=2)

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig5b-c1.pdf", width = 10,height = 4)
p0
dev.off()

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig5b-c2.pdf", width = 8,height = 4)
p1
dev.off()

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig5d.pdf", width = 16,height = 4)
p2
dev.off()

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig5d2.pdf", width = 8,height = 8)
p3
dev.off()

# Figure S7 
pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/figS7.pdf", width = 16,height = 4)
cowplot::plot_grid(plotlist = list(GSE115978_plotlist$VlnPlot_markers,GSE72056_plotlist$VlnPlot_markers), ncol=2)
dev.off()
