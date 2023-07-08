rm(list=ls())
setwd("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival")

pvalue <- list()
pvalue[["primary"]] <- openxlsx::read.xlsx("geneExprCorrelation/PRGs_74.p_value.xlsx", sheet = 1)
pvalue[["metastatic"]] <- openxlsx::read.xlsx("geneExprCorrelation/PRGs_74.p_value.xlsx", sheet = 2)

cor_list <- list()
cor_list[["primary"]] <- openxlsx::read.xlsx("geneExprCorrelation/PRGs_74.exprCor.xlsx", sheet = 1)
cor_list[["metastatic"]] <- openxlsx::read.xlsx("geneExprCorrelation/PRGs_74.exprCor.xlsx", sheet = 2)

barplot_data_list <- list()
for(tumor_type in c("primary","metastatic")){
  barplot_data <- data.frame()
  pvalue_tmp <- pvalue[[tumor_type]]
  for(gene in colnames(pvalue_tmp)[-1]){
    tmp <- data.frame(cor = cor_list[[tumor_type]][[gene]],p = pvalue_tmp[[gene]])
    tmp$group <- ifelse(tmp$cor >= 0.3, "Positive", ifelse(tmp$cor <= -0.3, "Negative", "None"))
    tmp$group[tmp$p >= 0.05] <- "None"
    tmp2 <- data.frame(gene = rep(gene,2), correlation = c("Positive","Negative"), value = c(sum(tmp$group == "Positive")-1,sum(tmp$group == "Negative")))
    barplot_data <- rbind(barplot_data,tmp2)
    barplot_data
  }
  barplot_data_list[[tumor_type]] <- barplot_data
}

barplot_data_list$primary$tumor_type <- "Primary"
barplot_data_list$metastatic$tumor_type <- "Metastatic"
mydata <- rbind(barplot_data_list$primary,barplot_data_list$metastatic)
mydata$value[mydata$correlation == "Negative"] <- 0-mydata$value[mydata$correlation == "Negative"]
mydata$tumor_type <- factor(mydata$tumor_type, levels = c("Primary", "Metastatic"))


# 转移原发分面
pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/figS1-barplot.pdf",width = 10,height = 9)
ggplot(data = mydata) +
  facet_wrap(~ tumor_type) +
  geom_bar(aes(x = gene, y = value, fill = correlation), width = 0.5, stat="identity")+
  labs(x = "", y = "") +
  coord_flip()+
  scale_fill_manual(values = c("#0072B5FF","#BC3C29FF"))+
  theme_bw()
dev.off()

