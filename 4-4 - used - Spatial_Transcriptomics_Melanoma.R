# 


rm(list=ls())
getwd()
# [1] "d:/Work/Projects/BMC_Revised/Results/04.singleCell_ST"
library(BayesSpace)
library(SingleCellExperiment)
library(ggplot2)
library(tidyr)
library(Matrix)
library(GSVA)

load("d:/Projects/Melanoma_Pyroptosis_new/8.spatial_transcriptome/RData/melanoma_and_spatialEnhancedMelanoma.RData")
# melanoma <- getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")
Gene_lists <- list(PScore = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])
Gene_lists$PScore[Gene_lists$PScore == "NLRP6"] <- "PAN3"


# ============================================================================= #
# 1.compute score ####
# ============================================================================= #
data <- logcounts(melanoma) %>% as.matrix()
data_E <- GSVA::gsva(data,Gene_lists,method="ssgsea")
data_ES <- t(data_E) %>% data.frame()
head(data_ES)
#           PScore
# 7x15  0.24491697
# 7x16  0.06023407
# 7x17  0.14672735
# 7x18 -0.02178160
# 8x13  0.31839771
# 8x14  0.58264378
save(data_ES, file = "d:/Work/Projects/BMC_Revised/Results/04.singleCell_ST/data_ES_spatial.RData")

data_ES_scale <- scale(data_ES) 
head(data_ES_scale)
#          PScore
# 7x15  0.7739341
# 7x16 -0.4269888
# 7x17  0.1354442
# 7x18 -0.9603057
# 8x13  1.2517516
# 8x14  2.9700438


# ============================================================================= #
# 2. volin plot ####
# ============================================================================= #
library(ggpubr)

colData <- melanoma@colData %>% data.frame()
colData <- cbind(colData,data_ES)
colData$defineTypes <- ifelse(colData$cluster.init==1, "Macrophage", 
                              ifelse(colData$cluster.init==2, "Stromal", 
                                     ifelse(colData$cluster.init==3, "Melanoma", "T_B")))
colData$defineTypes <- factor(colData$defineTypes,levels = c("Melanoma","Stromal","T_B","Macrophage"))
table(colData$defineTypes)

head(colData)
#      row col sizeFactor cluster.init spatial.cluster      PScore defineTypes
# 7x15   7  15  0.7955876            1               1  0.24491697  Macrophage
# 7x16   7  16  0.3073036            1               1  0.06023407  Macrophage
# 7x17   7  17  0.3312472            2               2  0.14672735     Stromal
# 7x18   7  18  0.4207466            3               2 -0.02178160    Melanoma
# 8x13   8  13  0.2554533            1               1  0.31839771  Macrophage
# 8x14   8  14  1.4734389            1               1  0.58264378  Macrophage


compare <- list(c("Stromal","Melanoma"),c("Melanoma","T_B"),c("Macrophage","Melanoma"),
                c("Stromal","T_B"),c("Macrophage","Stromal"),c("Macrophage","T_B"))

violinPlot <-
  ggplot(data = colData,aes(x=defineTypes,y=PScore,fill=defineTypes))+
  geom_violin(width=0.8,trim = F)+
  geom_boxplot(width = 0.2)+
  scale_x_discrete(limits=c("Melanoma","Stromal","T_B","Macrophage"),labels=c("Melanoma(95)","Stromal(121)","T/B(36)","Macrophage(41)"))+
  # ggsci::scale_fill_nejm() +
  # ggsci::scale_color_nejm() +
  scale_fill_manual(limits=c("Melanoma","Stromal","T_B","Macrophage"),values = c("#808080","#A000A0","#00FFFF","gold"),name="Response",guide= "none")+
  stat_compare_means(comparisons = compare,method = "wilcox.test")+
  xlab("") + ylab("Pyroptosis score")+ labs(title = "Spatial Transcription")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(color="black",size=12,angle = 0,hjust = 0.5,vjust = 1),
          axis.text.y = element_text(color="black",size=10),
          axis.title.y = element_text(color="black",size=12),
          axis.ticks.x = element_blank(),legend.position = "top",legend.direction = "horizontal",
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          strip.background = element_blank())
violinPlot

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig5g.pdf",width = 6,height = 5)
violinPlot
dev.off()


# ============================================================================= #
# 3.PScore in ST ####
# ============================================================================= #

logcounts_initial <- logcounts(melanoma) %>% as.matrix() %>% as.data.frame()

# merge
data_ES_t <- data_ES %>% t()
logcounts_add <- rbind(logcounts_initial,data_ES_t) %>% as.matrix()
logcounts_add_dgc <- as(logcounts_add, "dgCMatrix")

# new SingleCellExperiment
melanoma_PScore_sce <- SingleCellExperiment(assays=list(logcounts=logcounts_add_dgc), 
                                            colData=melanoma@colData)
set.seed(102)
melanoma_PScore_sce <- spatialPreprocess(melanoma_PScore_sce, platform="ST", n.PCs=7, n.HVGs=2000, log.normalize=F)

# cluster
melanoma_PScore_sce <- qTune(melanoma_PScore_sce, qs=seq(2, 10), platform="ST", d=7)
qPlot(melanoma_PScore_sce)
set.seed(149)
melanoma_PScore_sce <- spatialCluster(melanoma_PScore_sce, q=4, platform="ST", d=7,
                       init.method="mclust", model="t", gamma=2,
                       nrep=10000, burn.in=100,
                       save.chain=TRUE)


# spatialEnhance
melanoma_PScore_sce.enhanced <- spatialEnhance(melanoma_PScore_sce, q=4, platform="ST", d=7,
                                       model="t", gamma=2,
                                       jitter_prior=0.3, jitter_scale=3.5,
                                       nrep=100000, burn.in=100,
                                       save.chain=TRUE)

# this step is 20-30 mins, save data
save(melanoma_PScore_sce.enhanced,file = "melanoma_PScore_sce.enhanced.RData")

load("RData/melanoma_PScore_sce.enhanced.RData")
markers <- "PScore"
sce.enhanced.PScore <- enhanceFeatures(melanoma_PScore_sce.enhanced, melanoma_PScore_sce,
                                           feature_names=markers,
                                           nrounds=0)

logcounts(sce.enhanced.PScore)[markers, 1:5]
rowData(sce.enhanced.PScore)[markers, ]

# picture
enhanced.plots_1 <- featurePlot(sce.enhanced.PScore, markers)

enhanced.plots_2 <-featurePlot(sce.enhanced.PScore, markers)+
  scale_fill_gradient2(low = "#56B1F7",   mid="white",  high = "purple",
                       midpoint = 0.3,
                       space = "Lab", na.value = "#FFFFFF00", 
                       guide = "colourbar", aesthetics = "fill" )

clusterplot <- clusterPlot(melanoma_PScore_sce.enhanced, color="black")+
  scale_fill_manual(limits=c("1","2","3","4"),values=c("gold","#A000A0","#808080","#00FFFF"),label = c("1 Macrophage", "2 Stromal", "3 Melanoma", "4 T/B Cells"))


pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig5f.pdf", width = 14,height = 7)
patchwork::wrap_plots(list(clusterplot,enhanced.plots_2,enhanced.plots_1), ncol=3)
dev.off()


# ============================================================================= #
# 4. genes in ST ####
# ============================================================================= #
load("d:/Projects/Melanoma_Pyroptosis_new/8.spatial_transcriptome/RData/melanoma_and_spatialEnhancedMelanoma.RData")

Gene_lists <- list(PScore = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])
Gene_lists$PScore[Gene_lists$PScore == "NLRP6"] <- "PAN3"

data <- logcounts(melanoma) %>% as.matrix()
Gene_lists$PScore[Gene_lists$PScore %in% rownames(data)==FALSE]

# enhanceFeatures
markers <- Gene_lists$PScore
melanoma.enhanced_genes <- enhanceFeatures(melanoma.enhanced, melanoma,
                                     feature_names=markers,
                                     nrounds=0)

logcounts(melanoma.enhanced_genes)[markers, 1:5]
marker_ordered <- sort(markers)

# plot
# default color
enhanced.plots_1 <- purrr::map(marker_ordered, function(x){
  featurePlot(melanoma.enhanced_genes, x)})

p1 <- patchwork::wrap_plots(enhanced.plots_1, ncol=5)
p1

# change color
enhanced.plots_2 <- purrr::map(marker_ordered, function(x){
  featurePlot(melanoma.enhanced_genes, x)+
    scale_fill_gradient2(low = "#56B1F7",   mid="white",  high = "purple",
                         midpoint = c(max(logcounts(melanoma.enhanced_genes)[x,])/2+min(logcounts(melanoma.enhanced_genes)[x,])/2),
                        space = "Lab", na.value = "#FFFFFF00",   
                        guide = "colourbar", aesthetics = "fill" )})
p2 <- patchwork::wrap_plots(enhanced.plots_2, ncol=5)
p2


pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/figS8.pdf", width = 14,height = 19.8)
p2
dev.off()


