# 
# protein
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/05.ICB survival and response/"))
library(dplyr)
library(ggplot2)
library(ggpubr)
load((paste0(path, "Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")))
mela_exper <- read.table("d:/Projects/Melanoma_Pyroptosis/Genelist/melanoma_mRNA.txt")  #FPKM
ClinicalData <- read.delim("d:/Work/Data/public_protein/Melanoma_Immunotherapy_Protein_2021/Melanoma_protein_2021_ClinicalData.txt")
Protein <- read.delim("d:/Work/Data/public_protein/Melanoma_Immunotherapy_Protein_2021/Melanoma_protein_2021_FilterData60%.txt")


# ============================================================================ #
# 0.data prepare ####
# ============================================================================ #
colnames(Protein)
Protein <- Protein[Protein$Gene.names!="",]
Protein$Symbol <- data.frame(do.call(rbind,strsplit(Protein$Gene.names,"\\;")))$X1  # Protein$Gene.names[245] "WASH3P;WASH2P;WASH1"

Protein_f <- Protein[,c(ncol(Protein),1:174)]
table(duplicated(Protein_f$Symbol))
Protein_m <- as.data.frame(t(sapply(split(Protein_f,Protein_f$Symbol),
                                    function(x) apply(x[,2:ncol(x)],2,median))))
Protein_m[1:4,1:4]
data <- Protein_m %>% as.matrix()
data_t <- data %>% t() %>% data.frame()
data_t$sample <- rownames(data_t)


genelist <- list(PRGs_prot31_8 = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"])
table(rownames(Protein_m) %in% mela_exper$symbol)
# FALSE  TRUE 
#   266  4562 
proteins <- rownames(Protein_m)[rownames(Protein_m) %in% mela_exper$symbol ==FALSE]
table(genelist$PRGs_prot31_8 %in% rownames(Protein_m))
# FALSE  TRUE 
# 23     8 
genelist$PRGs_prot[genelist$PRGs_prot %in% rownames(Protein_m)]
# 8/31
# [1] "APIP"   "CASP1"  "CASP3"  "CASP4"  "CASP8"  "CHMP2B" "CHMP5"  "GSDMD" 
genelist$PRGs_prot <- c("APIP","CASP1", "CASP3", "CASP4", "CASP8", "CHMP2B", "CHMP5", "GSDMD")


# ============================================================================ #
# 1.compute score ###
# ============================================================================ #
data_ES <- GSVA::gsva(data,genelist,method="ssgsea")
data_ES <- t(data_ES) %>% data.frame()
data_ES$sample <- rownames(data_ES)

# save(data_ES,Protein_m,file = "Protein_PyropScore.Rdata")
# load("Protein_PyropScore.Rdata")

# merge clinical data
Clinical_PyropScore <- merge(data_ES,ClinicalData,by.x="sample",by.y="Sample.name")
Clinical_PyropScore$OS_Status <- ifelse(Clinical_PyropScore$Death=="N",0,1) 


# ============================================================================ #
# 2.boxplot ####
# ============================================================================ #
PScore_method = "PRGs_prot"

Clinical_PyropScore_F <- Clinical_PyropScore[Clinical_PyropScore$Immunotherapy.group%in%c("a-PD1"),]
Clinical_PyropScore_F$PyropScore_Scale <- scale(Clinical_PyropScore_F[[PScore_method]])

Clinical_PyropScore_F$Response2 <- ifelse(Clinical_PyropScore_F$Response%in%c("CR","PR"),"R","NR")
compare=list(c("R","NR"))

boxplot <- 
  ggplot(data = Clinical_PyropScore_F,aes(x=Response2,y=PyropScore_Scale,fill=Response2))+
  geom_boxplot(width = 0.6)+
  scale_x_discrete(limits=c("R","NR"),labels=paste(c("PRCR","SDPD"),"(",table(Clinical_PyropScore_F$Response2)[c(2,1)],")"))+
  scale_fill_manual(limits=c("R","NR"),values = c("#BC3C29FF","#0072B5FF"),name="Response")+
  stat_compare_means(comparisons = compare,method = "wilcox.test")+
  xlab("") + ylab("Pyroptosis Score")+ labs(title = PScore_method)+
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
boxplot

# ============================================================================ #
# 3. Survival and chisq####
# ============================================================================ #
library(survival)
library(survminer)

melcli_group <- Clinical_PyropScore_F
melcli_group$PyropScore <- as.numeric(melcli_group[[PScore_method]])
melcli_group$OS_use <- as.numeric(melcli_group$OS_Status)
melcli_group$OS_Time_use <- as.numeric(melcli_group$Overall.survival..months.)
melcli_group <- melcli_group[which(melcli_group$OS_use!="" & melcli_group$OS_Time_use!="" ),] 

res.cut <- surv_cutpoint(melcli_group, time = "OS_Time_use", event = "OS_use",variables = c("PyropScore"))
res.cat <- surv_categorize(res.cut)
melcli_group$Group <- ifelse(res.cat$PyropScore == "high","High_PS","Low_PS")
melcli_group$Group <- factor(melcli_group$Group,levels = c("Low_PS","High_PS"))
table(melcli_group$Group)

diff=survdiff(Surv(OS_Time_use,OS_use) ~ Group,data = melcli_group)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue=format(pValue, scientific = TRUE)
pValue

fit <- survfit(Surv(OS_Time_use,OS_use) ~ Group, data = melcli_group)
  
survival_plot <- 
    ggsurvplot(fit,
               risk.table=TRUE, 
               conf.int=F, 
               legend.labs=c("Low_PS","High_PS"), 
               legend.title="",
               xlab="Time(months)",
               ylab="Overall Survival",
               title=paste0("Proteomics data\n",PScore_method),
               palette = c("#0072B5FF","#BC3C29FF"),#颜色设置
               pval=paste0("p=",pValue),
               pval.method=T)
survival_plot 


tmp <- table(melcli_group$Group,melcli_group$Response2)
if(min(tmp)>=5){
  chisq_p <- chisq.test(melcli_group$Group,melcli_group$Response2)
  pvalue <- format(signif(chisq_p$p.value,2), scientific = TRUE)
}else{
  fisher_p <- fisher.test(melcli_group$Group,melcli_group$Response2)
  pvalue <- format(signif(fisher_p$p.value,2), scientific = TRUE)
}
melcli_group$Response2 <- factor(melcli_group$Response2, levels = c("R","NR"))

barplot <- 
  ggplot(data=melcli_group, mapping=aes(x=Group,fill=Response2))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c("#BC3C29FF","#0072B5FF"))+ 
  geom_text(stat='count',aes(label=..count..)
            , color="white", size=3.5,position=position_fill(0.5))+
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
          strip.background = element_blank())+
  labs(title=paste0("Proteomics data\n",PScore_method), subtitle = paste0("p value = ",pvalue) , y="Pecentage(%)")
barplot


pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig6f-g.pdf",width = 10,height = 7)
cowplot::plot_grid(plotlist = list(boxplot,barplot),nrow=1)
dev.off()

pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/fig6e.pdf",onefile = F,width = 5,height =7)
survival_plot
dev.off()


