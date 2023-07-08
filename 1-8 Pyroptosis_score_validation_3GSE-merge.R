# pyroptosis_related_datasets, merged figure


rm(list=ls())
options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))

load("PScoreValidation/GSE57253_ggplot.RData")
load("PScoreValidation/GSE153494_ggplot.RData")
load("PScoreValidation/GSE192714_ggplot.RData")

pyroptosis_related_datasets <- rbind(GSE57253_ggplot,GSE153494_ggplot,GSE192714_ggplot)



p1 <- pyroptosis_related_datasets %>% 
        as_tibble() %>% 
        dplyr::select(PS:geo_accession) %>%
        mutate(geo_accession = factor(geo_accession, levels = c("GSE57253","GSE153494","GSE192714"))) %>%   # 对数据集和图例进行排列
        group_by(geo_accession) %>%
        ggplot(aes(x=group,y=PS,color=group))+
          facet_grid(~ geo_accession,scales = "free", space = "free") + 
          geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+
          ggbeeswarm::geom_beeswarm(alpha=0.5,size=0.5) + 
          theme_bw()
p1

data_text<-data.frame(label=c("Kruskal-Wallis, p = 5.7e-3", "Kruskal-Wallis, p = 1.4e-2", "Wilcoxon, p = 2.9e-2"),
                      geo_accession = factor(c("GSE57253","GSE153494","GSE192714")),
                      group = c("Healthy Controls","10m","Control"),
                      PS = c(1.25,1.25,1.25))

#                         label geo_accession            group   PS
# 1 Kruskal-Wallis, p = 5.7e-3      GSE57253 Healthy Controls  1.25
# 2 Kruskal-Wallis, p = 1.4e-2     GSE153494             10m   1.25
# 3       Wilcoxon, p = 2.9e-2     GSE192714          Control  1.25

p1 + geom_text(data=data_text,mapping = aes(x=group,y=PS,label=label),
               nudge_x=0.5,nudge_y=0.2,size = 3,colour = "black")

