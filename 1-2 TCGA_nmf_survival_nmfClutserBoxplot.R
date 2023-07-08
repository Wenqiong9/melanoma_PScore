# Aim: dividing TCGA melanoma patients into 2 cluster by NMF
# Aim: survival analysis
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
library(survival)
library(survminer)
path <- "d:/Work/Projects/BMC_Revised/"
setwd(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/"))
load("coxForest/melanoma_expr_list.PRGs_74.RData")
load("coxForest/Expr_Survival_list.PRG_74.RData")

get_nmf_res <- function(input=NULL, outfile=NULL, rank_space=2:6, seed = NULL, nrun = 100, thread=24, cluster_res=FALSE){
  library(NMF)
  # input: row-feature, col-samples
  nmf_res <- nmf(x=input, rank=rank_space, seed=seed, nrun = nrun, .options=paste0("p", thread))
  pdf(outfile)
  print(plot(nmf_res))
  dev.off()
  # We select values of k where the magnitude of the cophenetic correlation coefficient begins to fall
  # Ref: Metagenes and molecular pattern discovery using matrix factorization
  coph <- nmf_res$measures$cophenetic
  coph_diff <- coph[-length(coph)]-coph[-1]
  rank_opt <- which.max(coph_diff)+1
  if(cluster_res){
    opt_rank_name <- as.character(rank_opt)
    nmf_res_opt <- nmf_res$fit[[opt_rank_name]]
    feat_opt_idx <- extractFeatures(object=nmf_res_opt, method = "max", format = "list", nodups = TRUE)
    feat_opt_idx <- as.numeric(na.omit(unique(unlist(feat_opt_idx))))
    feat_opt <- rownames(input)[feat_opt_idx]
    nmf_cluster_opt <- nmf(input[feat_opt, ], rank = as.numeric(rank_opt), seed = seed)
    group <- predict(nmf_cluster_opt)
    data_group <- data.frame(Sample=names(group), Group=group)
    return(list(nmf_res=nmf_res, rank_opt_coph=rank_opt, feat_opt=feat_opt, group_opt=data_group))
  }else{
    return(list(nmf_res=nmf_res, rank_opt_coph=rank_opt))
  }
}



# ============================================================================ #
# 1.NMF -- metastatic patients ####
# ============================================================================ #
#### 1.1 NMF ####
GeneExpr <- melanoma_expr_list[["metastatic"]]
GeneExpr <- GeneExpr[,gsub("\\.", "-", substr(colnames(GeneExpr), 1, 12)) %in% rownames(Expr_Survival_list[["metastatic"]])]

# NMF
i <- 6666
nmf_metastatic <- get_nmf_res(GeneExpr,outfile = "nmf_survival/NMF_metastatic_RPGs74_cophenetic.pdf",
                   rank_space=2:7,nrun = 70, thread=24,cluster_res=TRUE,seed = NULL)
save(nmf_metastatic,file = "nmf_survival/NMF_metastatic_RPGs74_Res.Rdata")
group_metastatic = nmf_metastatic$group_opt
print(table(group_metastatic$Group))
#   1   2   3 
# 149  51 161 

pdf("nmf_survival/NMF_metastatic_RPGs74_Heatmap.pdf",width = 5,height = 7)
consensusmap(nmf_metastatic$nmf_res$consensus$`3`,labCol=NA, labRow=NA)
dev.off()


#### 1.2 survival analysis of 2 NMF cluster ####

# merge NMF cluster and clinic information 
load("nmf_survival/NMF_metastatic_RPGs74_Res.Rdata")
group_metastatic = nmf_metastatic$group_opt

group_metastatic$Sample <- gsub("\\.","-",substr(group_metastatic$Sample,1,12))
mel_cli <- Expr_Survival_list$metastatic
mel_cli$Sample <- rownames(mel_cli)
melcli_group <- merge(group_metastatic, mel_cli, by = "Sample") 
melcli_group$OS_year <- melcli_group$OS.time/365

# survival analysis
diff = survdiff(Surv(OS.time, OS) ~ Group,data = melcli_group)
pValue = 1-pchisq(diff$chisq,df=1)
pValue = signif(pValue,2)
pValue
pValue = format(pValue, scientific = TRUE)

# HR
cox.res=coxph(Surv(OS.time,OS)~ Group,data = melcli_group)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)
CI_group2_low <- signif(CI[1],2)
CI_group2_high <- signif(CI[3],2)
CI_group3_low <- signif(CI[2],2)
CI_group3_high <- signif(CI[4],2)

fit <- survfit(Surv(OS_year, OS) ~ Group, data = melcli_group)

ggsurvplot(fit, 
           data = melcli_group,
           conf.int=F,
           pval= paste0("p = ",pValue,
                       "\nHR (Cluster 2) = ",HR[1]," (",CI_group2_low," - ",CI_group2_high,")",
                       "\nHR (Cluster 3) = ",HR[2]," (",CI_group3_low," - ",CI_group3_high,")"),
           risk.table=TRUE,
           legend.labs=c("Cluster 1", "Cluster 2","Cluster 3"),
           legend.title="Cluster")


#### 1.3 computePScore ####
# 1-4 computePScore.R
load(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/PScoreModel_ROC/PyropScore_SKCM.Rdata"))
data_ES_zscore <- data_ES_zscore_List$metastatic


#### 1.4 boxplot ####
library(ggplot2)
library(ggpubr)

# load("nmf_survival/NMF_metastatic_RPGs74_Res.Rdata")
# group_metastatic = nmf_metastatic$group_opt
data_ES_zscore$Sample <- rownames(data_ES_zscore)
data_ES_zscore_merge <- merge(group_metastatic,data_ES_zscore,by="Sample")
data_ES_zscore_merge$Group <- factor(data_ES_zscore_merge$Group)

my_comparisons <- list(c("1","2"),c("2","3"),c("1","3"))
ggplot(data_ES_zscore_merge,aes(x=Group,y=.data[["PRGs_prot"]],color=Group))+
  geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+
  ggbeeswarm::geom_beeswarm(size=0.5)+
  labs(y="PScore")+
  ggpubr::stat_compare_means(comparisons=my_comparisons) +
  theme_classic()



# ============================================================================ #
# 2.NMF -- primary patients  ####
# ============================================================================ #
#### 2.1 NMF 3 clusters ####
GeneExpr <- melanoma_expr_list[["primary"]]
GeneExpr <- GeneExpr[,gsub("\\.", "-", substr(colnames(GeneExpr), 1, 12)) %in% rownames(Expr_Survival_list[["primary_OS"]])]

i = 6666
nmf_primary <- get_nmf_res(GeneExpr,outfile = "nmf_survival/NMF_primary_PRGs74_cophenetic.pdf",
                   rank_space=2:7,nrun = 70, thread=24,cluster_res=TRUE,seed = i)
save(nmf_primary,file = "nmf_survival/NMF_primary_RPGs74_Res.Rdata")

pdf("nmf_survival/NMF_primary_RPGs74_Heatmap.pdf",width = 5,height = 7)
consensusmap(nmf_primary$nmf_res$consensus$`3`,labCol=NA, labRow=NA)
dev.off()


#### 2.2 survival analysis of 2 NMF cluster ####

# merge NMF cluster and clinic information 
load("nmf_survival/NMF_primary_RPGs74_Res.Rdata")
group_primary = nmf_primary$group_opt
print(table(group_primary$Group))
#  1  2  3 
# 19 56 28 

group_primary$Sample <- gsub("\\.","-",substr(group_primary$Sample,1,12))
mel_cli <- Expr_Survival_list$primary
mel_cli$Sample <- rownames(mel_cli)
melcli_group <- merge(group_primary, mel_cli, by = "Sample") 

# survival 
diff=survdiff(Surv(OS.time, OS) ~ Group,data = melcli_group)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue
pValue=format(pValue, scientific = TRUE)

# HR
cox.res=coxph(Surv(OS.time,OS)~ Group,data = melcli_group)
HR <- signif(exp(coef(cox.res)),2)
CI <- signif(exp(confint(cox.res)),2)
CI_group2_low <- signif(CI[1],2)
CI_group2_high <- signif(CI[3],2)
CI_group3_low <- signif(CI[2],2)
CI_group3_high <- signif(CI[4],2)

melcli_group$OS_year <- melcli_group$OS.time/365
fit <- survfit(Surv(OS_year, OS) ~ Group, data = melcli_group)

ggsurvplot(fit, 
           data = melcli_group,
           conf.int=F,
           pval= paste0("p = ",pValue,
                        "\nHR (Cluster 2) = ",HR[1]," (",CI_group2_low," - ",CI_group2_high,")",
                        "\nHR (Cluster 3) = ",HR[2]," (",CI_group3_low," - ",CI_group3_high,")"),
           risk.table=TRUE,
           legend.labs=c("Cluster 1", "Cluster 2","Cluster 3"),
           legend.title="Cluster")


# #### 2.3 computePScore ####
# 1-4 computePScore.R
load(paste0(path, "Results/01.heatmap_NMF_coxForest_survival/PScoreModel_ROC/PyropScore_SKCM.Rdata"))
data_ES_zscore <- data_ES_zscore_List$primary

#### 2.4 boxplot ####
library(ggplot2)
library(ggpubr)

load("nmf_survival/NMF_primary_RPGs74_Res.Rdata")
group_primary = nmf_primary$group_opt
data_ES_zscore$Sample <- rownames(data_ES_zscore)
data_ES_zscore_merge <- merge(group_primary,data_ES_zscore,by="Sample")
data_ES_zscore_merge$Group <- factor(data_ES_zscore_merge$Group)
 
PScore_method = names(data_ES_zscore)[4]  # PRGs_prot

my_comparisons <- list(c("1","2"),c("2","3"),c("1","3"))

ggplot(data_ES_zscore_merge,aes(x=Group,y=.data[["PRGs_prot"]],color=Group))+
  geom_boxplot(width=0.6,fill=NA,outlier.shape = NA)+
  ggbeeswarm::geom_beeswarm(size=0.5)+
  labs(y="PScore")+
  ggpubr::stat_compare_means(comparisons=my_comparisons) +
  theme_classic()


