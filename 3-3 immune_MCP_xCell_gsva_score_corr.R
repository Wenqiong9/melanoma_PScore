

rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
my.cor.test <- function(...) {
  obj<-try(cor.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/03.immuneMechanism/"))
load("gset_all.add_PScore.Rdata")
load("score_res_mcp.RData")


# ============================================================================ #
# correlation ####
# ============================================================================ #
str(score_res)
gset_all$GEO_ID == names(score_res)
dataset_id <- gset_all$GEO_ID


mcp_cor_list <- list()
for(j in 1:7){
  PScore <- gset_all$PScore[[j]]
  PScore$ID <- rownames(PScore)
  PScore <- PScore[,c("ID","PScore")]
  
  Data <- score_res[[j]] %>% data.frame()
  Data <- merge(PScore, Data, by="ID")
  
  Cor <- data.frame(Gene=colnames(Data)[!colnames(Data) %in% c("ID","PScore")])
  Cor[,c("estimate.rho","pvalue")] <- t(apply(Data[Cor$Gene], 2, function(x)unlist(my.cor.test(as.numeric(x),as.numeric(Data$PScore),method="pearson"))))  # spearman
  Cor$FDR <- p.adjust(Cor$pvalue,method = "fdr")
  
  mcp_cor_list[[dataset_id[j]]] <- Cor
}
 
save(mcp_cor_list, file = "mcp_cor_list.RData")


# ============================================================================ #
# plot ####
# ============================================================================ #
library(ggplot2)

immuneScore_PScore_Corr <- mcp_cor_list %>% dplyr::bind_rows() %>% dplyr::mutate(Dataset=rep(gset_all$GEO_ID, 10))  # cell types of diff methods
immuneScore_PScore_Corr$Gene <- gsub("_MCPcounter","",immuneScore_PScore_Corr$Gene)
immuneScore_PScore_Corr$FDR_log10 <- -log10(immuneScore_PScore_Corr$FDR)
immuneScore_PScore_Corr$FDR_log10[immuneScore_PScore_Corr$FDR_log10 > -log(0.05)] <- 3

# plotting parameters 
immuneScore_PScore_Corr$class <- ifelse(immuneScore_PScore_Corr$estimate.rho>0,"Pos","Neg")
y_order <- t(sapply(split(immuneScore_PScore_Corr[,c("estimate.rho","FDR_log10")],immuneScore_PScore_Corr$Gene),function(x) 
  c(nrow(x[x$estimate.rho > 0 & x$FDR_log10 > -log10(0.05),]),nrow(x[x$estimate.rho < 0 & x$FDR_log10 > -log10(0.05),])))) %>% 
  data.frame() %>% 
  magrittr::set_colnames(c("Pos","Neg"))
y_order <- y_order[order(y_order$Pos-y_order$Neg),]

estimate_range <- max(abs(range(immuneScore_PScore_Corr$estimate.rho,na.rm = T)))

pdf("immune_mechanism_mcp.pdf",width = 5, height = 3.5)
ggplot(immuneScore_PScore_Corr,aes(x=Gene,y=Dataset,color=estimate.rho,fill=estimate.rho))+
  geom_point(aes(size=FDR_log10),shape=23)+
  scale_fill_gradientn(limit=c(-estimate_range,estimate_range),colors= colorRampPalette(c("#0072B5FF","white","red"),space="rgb")(100),breaks=seq(-1,1,length.out = 5),labels=c("-1","-0.5","0","0.5","1"),name="Rs")+
  scale_color_gradientn(limit=c(-estimate_range,estimate_range),colors= colorRampPalette(c("#0072B5FF","white","red"),space="rgb")(100),breaks=seq(-1,1,length.out = 5),labels=c("-1","-0.5","0","0.5","1"),name="Rs")+
  scale_size_continuous(limit=c(0,3),range = c(0.05, 3),breaks=c(-log10(0.05),-log10(0.01),-log10(0.001)),labels=c("0.05","0.01","1e-3"),name="FDR")+
  geom_point(data=immuneScore_PScore_Corr[immuneScore_PScore_Corr$FDR_log10 > -log10(0.05),],aes(y=Dataset,x=Gene,size=FDR_log10),shape=23,color="black") +
  coord_fixed()+
  scale_y_discrete(limit=rev(names(table(immuneScore_PScore_Corr$Dataset))))+
  scale_x_discrete(limit=rownames(y_order),label=gsub("\\."," ",rownames(y_order)))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(size=8,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),legend.position = "bottom",
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()
