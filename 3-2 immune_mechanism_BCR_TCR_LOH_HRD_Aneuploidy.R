# TCR、BCR、 Aneuploidy  LOH_frac_altered  HRD -- from panimmune
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
# TCR_BCR_LOH_Aneuploidy_HRD ####
# ============================================================================ #
load("gset_all.add_PScore.Rdata")

TCR_BCR_Richness <- openxlsx::read.xlsx("d:/Work/Data/public_immune/TCGA - The Immune Landscape of Cancer/NIHMS958212-supplement-2.xlsx")
TCR_BCR_Richness <- TCR_BCR_Richness[TCR_BCR_Richness$TCGA.Study == "SKCM",]
TCR_Richness <- TCR_BCR_Richness[,c("TCR.Richness","TCGA.Participant.Barcode")]
BCR_Richness <- TCR_BCR_Richness[,c("BCR.Richness","TCGA.Participant.Barcode")]
TCR_Richness_log2 <- TCR_Richness
TCR_Richness_log2$TCR.Richness <- log2(TCR_Richness_log2$TCR.Richness +1)
BCR_Richness_log2 <- BCR_Richness
BCR_Richness_log2$BCR.Richness <- log2(BCR_Richness_log2$BCR.Richness +1)
(head(TCR_Richness));(head(BCR_Richness))

# "Silent.per.Mb", "Non.silent.per.Mb"
# mutationBurden <- read.delim("d:/Work/Data/public_immune/TCGA_panimmune/mutation-load_updated.txt")

# "LOH_n_seg", "LOH_frac_altered"
LOH <- readr::read_tsv("d:/Work/Data/public_immune/TCGA_panimmune/ABSOLUTE_scores.tsv")
LOH$TCGA.Participant.Barcode <- substr(LOH$...1,1,12)
LOH_n_seg <- LOH[,c("LOH_n_seg","TCGA.Participant.Barcode")]
LOH_n_seg_log2 <- LOH_n_seg
LOH_n_seg_log2$LOH_n_seg <- log2(LOH_n_seg$LOH_n_seg+1)
(head(LOH_n_seg))

Aneuploidy <- TCR_BCR_Richness[,c( "Aneuploidy.Score","TCGA.Participant.Barcode")]
HRD <- TCR_BCR_Richness[,c( "Homologous.Recombination.Defects","TCGA.Participant.Barcode")]    # TCGA.HRD_withSampleID.txt


SKCM_BCR_TCR_LOH_HRD_list <- list("log2(TCR)" = TCR_Richness_log2,
                                  "log2(BCR)" = BCR_Richness_log2,
                                  "log2(LOH_n_seg)"= LOH_n_seg_log2 ,
                                  Aneuploidy=Aneuploidy,
                                  HRD=HRD)

for(i in 1:length(SKCM_BCR_TCR_LOH_HRD_list)){
  temp <- SKCM_BCR_TCR_LOH_HRD_list[[i]]
  colnames(temp) <- c("Feature","barcode") 
  SKCM_BCR_TCR_LOH_HRD_list[[i]] <- temp
}



# ============================================================================ #
# picture #### 
# ============================================================================ #
data_ES_zscore <- gset_all$PScore[[2]]
data_ES_zscore$id <- gsub("\\.","-",substr(rownames(data_ES_zscore),1,12))
data_ES_zscore <- data_ES_zscore[substr(rownames(data_ES_zscore),14,15)<10,]
data_ES_zscore <- data_ES_zscore[!duplicated(data_ES_zscore$id),]


CorrPlot_all <- list()
for (i in names(SKCM_BCR_TCR_LOH_HRD_list)) {
  print(paste("Running for",i))
  
  sub <- SKCM_BCR_TCR_LOH_HRD_list[[i]]
  sub <- merge(sub,data_ES_zscore,by.x="barcode",by.y="id")
  CorrPlot_all[[i]] <- ggplot(sub, aes(x=PScore_zscore, y=Feature))+
    geom_point()+
    xlab("Pyroptosis Score")+
    ylab(i)+
    stat_smooth(method = "lm",se=T)+stat_cor(data=sub, method = "pearson")+
    theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())
}

pdf("immune_mechanism_TCR_BCR_LOH_Aneuploidy_HRD.pdf",width = 15,height = 3)
cowplot::plot_grid(plotlist = CorrPlot_all,ncol = 5)
dev.off()


