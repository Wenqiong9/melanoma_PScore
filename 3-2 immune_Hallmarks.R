# 
# Figure 3 hallmark
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
# 1.cancer hallmark####
# ============================================================================ #
load("gset_all.add_PScore.Rdata")
Hm <- fgsea::gmtPathways("d:/Work/Data/public_immune/Hallmark/h.all.v6.1.symbols.gmt")

for(i in 1:7){
  print(gset_all$GEO_ID[[i]])
  tmp <- gset_all$PScore[[i]]
  tmp$geo_accession <- rownames(tmp)
  gset_all$PScore[[i]] <- tmp
}

# hallmark score
gset_all <- gset_all %>% dplyr::mutate(HallmarkScore=purrr::map2(.x=gset_data_matrix_symbol,.y=PScore,function(.x,.y){
  try(
    { Data <- .x %>% as.matrix()
    Data_ES <- GSVA::gsva(Data,Hm)
    Data_ES <- Data_ES %>% t() %>% as.data.frame()
    Data_ES$geo_accession <- rownames(Data_ES)
    Data_ES <- merge(Data_ES,.y[c("PScore","geo_accession")],by="geo_accession")
    return(Data_ES)
    },silent = T)
}))

# cor
hallmak_Pyro_Corr <- purrr::map(.x=gset_all$HallmarkScore,function(.x){
  Data <- .x %>% data.frame()
  Cor <- data.frame(Gene=colnames(Data)[!colnames(Data)%in%c("geo_accession","PScore")])
  # pearson or spearman?
  Cor[,c("estimate.rho","pvalue")] <- t(apply(Data[Cor$Gene], 2, function(x)unlist(my.cor.test(as.numeric(x),as.numeric(Data$PScore),method="pearson"))))
  Cor$FDR <- p.adjust(Cor$pvalue,method = "fdr")
  return(Cor)
})

# save
Hallmark_PS_Corr <- hallmak_Pyro_Corr %>% dplyr::bind_rows() %>% dplyr::mutate(Dataset=rep(gset_all$GEO_ID,each=50))
write.csv(Hallmark_PS_Corr,file = "Hallmark_PS_Corr.csv",quote = F)


# ============================================================================ #
# 2.picture####
# ============================================================================ #
Hallmark_PS_Corr <- read.csv("Hallmark_PS_Corr.csv")

Hallmark_PS_Corr$FDR_log10 <- -log10(Hallmark_PS_Corr$FDR)
Hallmark_PS_Corr$FDR_log10[Hallmark_PS_Corr$FDR_log10 > -log(0.05)] <- -log10(1e-10)

Hallmark_PS_Corr$class <- ifelse(Hallmark_PS_Corr$estimate.rho>0,"Pos","Neg")
y_order <- t(sapply(split(Hallmark_PS_Corr[,c("estimate.rho","FDR_log10")],Hallmark_PS_Corr$Gene),function(x) 
  c(nrow(x[x$estimate.rho > 0 & x$FDR_log10 > -log10(0.05),]),nrow(x[x$estimate.rho < 0 & x$FDR_log10 > -log10(0.05),])))) %>% 
  data.frame() %>% 
  set_colnames(c("Pos","Neg"))
y_order <- y_order[order(y_order$Pos-y_order$Neg),]

estimate_range <- max(abs(range(Hallmark_PS_Corr$estimate.rho,na.rm = T)))


pdf("immune_mechanism_Hallmark.pdf",width=12,height = 4) 
ggplot(Hallmark_PS_Corr,aes(y=Dataset,x=Gene,color=estimate.rho,fill=estimate.rho))+
  geom_point(aes(size=FDR_log10),shape=23)+
  scale_fill_gradientn(limit=c(-estimate_range,estimate_range),colors= colorRampPalette(c("#0072B5FF","white","red"),space="rgb")(100),breaks=seq(-1,1,length.out = 5),labels=c("-1","-0.5","0","0.5","1"),name="Rs")+
  scale_color_gradientn(limit=c(-estimate_range,estimate_range),colors= colorRampPalette(c("#0072B5FF","white","red"),space="rgb")(100),breaks=seq(-1,1,length.out = 5),labels=c("-1","-0.5","0","0.5","1"),name="Rs")+
  scale_size_continuous(limit=c(0,10.1),range = c(0.1, 5),breaks=c(-log10(0.05),5,10),labels=c("0.05","1e-5","<1e-10"),name="FDR")+
  coord_fixed()+
  scale_y_discrete(limit=rev(names(table(Hallmark_PS_Corr$Dataset))))+
  scale_x_discrete(limit=rownames(y_order),label=Hmisc::capitalize(tolower(gsub("\\_"," ", gsub("HALLMARK_","",rownames(y_order))))))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(size=8,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),legend.position = "bottom",
        legend.key = element_rect(fill="white",colour = "black"))+
  geom_point(data=Hallmark_PS_Corr[Hallmark_PS_Corr$FDR_log10 > -log10(0.05),],aes(y=Dataset,x=Gene,size=FDR_log10),shape=23,color="black")
dev.off()

