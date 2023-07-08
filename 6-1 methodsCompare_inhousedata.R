# inhouse data
# 
# ============================================================================ #


rm(list=ls())
options(stringsAsFactors = F)
# path <- "d:/Work/Projects/BMC_Revised/"
# setwd(paste0(path, "Results/05.ICB survival and response/"))
library(dplyr)
library(ggplot2)
load("d:/Work/Projects/BMC_Revised/Results/01.heatmap_NMF_coxForest_survival/coxForest/uniCoxSig_list.RData")
load("d:/Work/Projects/BMC_Revised/Data/Gene_list/PRGs_74.RData")
inhouse_expr <- openxlsx::read.xlsx("d:/Work/Data/inhouse_data(xiangya)/Supplementary tables_HEYi_Ferroptosis.xlsx",sheet=9,startRow = 2,rowNames = TRUE)
inhouse_meta <- openxlsx::read.xlsx("d:/Work/Data/inhouse_data(xiangya)/Supplementary tables_HEYi_Ferroptosis.xlsx",sheet=8,startRow = 2)



# ============================================================================ #
# 1.data prepare -- inhouse data####
# ============================================================================ #

genelist <- list(PRGs_74 = PRGs_74$Symbol.ID,
                 PRGs_prot = uniCoxSig_list$metastatic$gene[uniCoxSig_list$metastatic$Risk_group == "Protect"],
                 Wang_Ann.Transl.Med_2022 = c("AIM2", "GSDMC", "GSDMD", "IL18", "NLRP6", "PRKACA"))

genelist$PRGs_74[genelist$PRGs_74 %in% rownames(inhouse_expr)==FALSE]
rownames(inhouse_expr)[rownames(inhouse_expr) == "PJVK"] <- "DFNB59"
rownames(inhouse_expr)[rownames(inhouse_expr)=="GSDME"] <- "DFNA5"
rownames(inhouse_expr)[rownames(inhouse_expr)=="LAMTOR3"] <- "MAPKSP1"
rownames(inhouse_expr)[rownames(inhouse_expr)=="LAMTOR4"] <- "C7orf59"
rownames(inhouse_expr)[rownames(inhouse_expr)=="LAMTOR5"] <- "HBXIP"
rownames(inhouse_expr)[rownames(inhouse_expr)=="IST1"] <- "KIAA0174"
rownames(inhouse_expr)[rownames(inhouse_expr)=="MAP3K20"] <- "ZAK"
rownames(inhouse_expr)[rownames(inhouse_expr)=="CHMP3"] <- "VPS24"
genelist$Wang_Ann.Transl.Med_2022[genelist$Wang_Ann.Transl.Med_2022 %in% rownames(inhouse_expr)==FALSE]

# ssgsea
expr_Data <- inhouse_expr %>% as.matrix()
data_ES <- GSVA::gsva(expr_Data,genelist,method="ssgsea")
data_ES <- data_ES %>% t() %>% data.frame()
head(data_ES)


# coefficient
inhouse_expr <- as.data.frame(t(inhouse_expr))

# Zhu_Stem.Cells.Int_2023: Score=(-0.166*CCL8)-0.156*FCGR2A+0.047*GBP2-0.327*GRIPAP1-0.207*HAPLN3+0.145*HPDL-0.022*IFITM1-1.01*TRIM34.
inhouse_expr$Zhu_Stem.Cells.Int_2023 <- (-0.166 * inhouse_expr$CCL8) - 0.156 * inhouse_expr$FCGR2A +
    0.047 * inhouse_expr$GBP2 - 0.327 * inhouse_expr$GRIPAP1 -
    0.207 * inhouse_expr$HAPLN3 + 0.145 * inhouse_expr$HPDL -  # HPDL: 4-HPPD-L, GLOXD1, NEDSWMA, SPG83
    0.022 * inhouse_expr$IFITM1 - 1.01 * inhouse_expr$TRIM34   # TRIM34: IFP1, RNF21
# Li_Int.J.Gen.Med_2022: Score=(-0.119*AIM2)+(-0.487*NLRP6)+(-0.374*IL18)+(0.230*GSDMA)+(0.383*GSDMC)
inhouse_expr$Li_Int.J.Gen.Med_2022 <- (-0.119 * inhouse_expr$AIM2) + (-0.487 * inhouse_expr$NLRP6)+
  (-0.374 * inhouse_expr$IL18) + (0.230 * inhouse_expr$GSDMA)+
  (0.383 * inhouse_expr$GSDMC)
# Xu_Front.Med_2022: Score=0.003*EMP3-0.065*TLR1-0.012*IFNGR2-0.288*IL15-0.057*CCL8-0.633*NLRP6-0.329*CCL25-0.024*RTP4
inhouse_expr$Xu_Front.Med_2022 <- (0.003 * inhouse_expr$EMP3) - 0.065 * inhouse_expr$TLR1 - 
  0.012 * inhouse_expr$IFNGR2 - 0.288 * inhouse_expr$IL15 - 
  0.057 * inhouse_expr$CCL8 - 0.633 * inhouse_expr$NLRP6 - 
  0.329 * inhouse_expr$CCL25 - 0.024 * inhouse_expr$RTP4
# Niu_Math.Biosci.Eng_2022: Score=0.038*GSDMA+0.31*GSDMC-0.028*AIM2-0.437*NOD2
inhouse_expr$Niu_Math.Biosci.Eng_2022 <- 0.038 * inhouse_expr$GSDMA + 0.31 * inhouse_expr$GSDMC - 
  0.028 * inhouse_expr$AIM2 - 0.437 * inhouse_expr$NOD2
# Wu_PeerJ_2021: Score=0.2758*GSDMC-0.0699*GZMA-0.0526*AIM2-0.1766*PDL1
inhouse_expr$Wu_PeerJ_2021 <- 0.2758 * inhouse_expr$GSDMC - 0.0699 * inhouse_expr$GZMA - 
  0.0526 * inhouse_expr$AIM2 - 0.1766 * inhouse_expr$CD274
# Ju_Front.Oncol_2021: Score=-0.006861*GSDMD+0.0003969*GSDME-0.001943*CASP4+0.0079361*GSDMC-0.022123*NLRC4-0.009636*APIP-0.003569*AIM2-0.00106*CASP3-0.000169*IL18
inhouse_expr$Ju_Front.Oncol_2021 <- -0.006861 * inhouse_expr$GSDMD + 0.0003969 * inhouse_expr$DFNA5 - 
  0.001943 * inhouse_expr$CASP4 + 0.0079361 * inhouse_expr$GSDMC - 
  0.022123 * inhouse_expr$NLRC4 - 0.009636 * inhouse_expr$APIP - 
  0.003569 * inhouse_expr$AIM2 - 0.00106 * inhouse_expr$CASP3 - 
  0.000169 * inhouse_expr$IL18
# Shi_Medicine_2022: Score=-0.0452*BST2-0.1636*GBP5-0.0531*AIM2
inhouse_expr$Shi_Medicine_2022 <-  -0.0452 * inhouse_expr$BST2 - 0.1636 * inhouse_expr$GBP5 - 0.0531 * inhouse_expr$AIM2
# Wu_Cancer.Med_2022: Score=0.139*CASP5+0.240*NLRP6+1.388*NLRP7+0.112*PYCARD
inhouse_expr$Wu_Cancer.Med_2022 <- 0.139 * inhouse_expr$CASP5 + 0.240 * inhouse_expr$NLRP6 + 
  1.388 * inhouse_expr$NLRP7 + 0.112 * inhouse_expr$PYCARD
# Wang_J.Investig.Dermatol_2022: Score = mean(IRF9, STAT2)
inhouse_expr$Wang_J.Investig.Dermatol_2022 <- apply(inhouse_expr[,c("IRF9","STAT2")],1,mean)


head(data_ES)
head(inhouse_expr[,59386:59394])
data_ES$ID <- rownames(data_ES)
data_ES$ID == rownames(inhouse_expr)

multiScore <- cbind(data_ES,inhouse_expr[,59386:59394])
multiScore$ID == inhouse_meta$ID
inhouse_merge <- merge(multiScore,inhouse_meta,by.x="ID")


# ============================================================================ #
# 2.cox forest ####
# ============================================================================ #
library(survival)
library(ggplot2)

uniCox <- data.frame()
# j = colnames(inhouse_merge)[3]
# j = colnames(inhouse_merge)[4]
# j = colnames(inhouse_merge)[11]
for(j in colnames(inhouse_merge)[2:13]){
  cox <- coxph(Surv(PFST.Month.,PFS) ~ inhouse_merge[,j], data = inhouse_merge)
  coxSummary <- summary(cox)
  C_index = signif(coxSummary$concordance[["C"]],3)
  uniCox <- rbind(uniCox, data.frame(gene = j,
                                     HR = coxSummary$conf.int[,"exp(coef)"],
                                     HR.95L = coxSummary$conf.int[,"lower .95"],
                                     HR.95H = coxSummary$conf.int[,"upper .95"],
                                     P_value = coxSummary$coefficients[,"Pr(>|z|)"],
                                     Concordance = C_index,
                                     Concordance_se = paste0(C_index, "(se = ", signif(coxSummary$concordance[["se(C)"]],3),")")))
}
uniCox$Risk_group <- ifelse(uniCox$HR >= 1 & uniCox$P_value < 0.05, "Risk",
                            ifelse(uniCox$HR < 1 & uniCox$P_value < 0.05, "Protect","Not_sig"))

ggplot(uniCox, aes(HR, gene)) + 
  geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
  geom_point(shape=18, aes(col=Risk_group, size = P_value))  +
  geom_vline(aes(xintercept = 1),color="darkgrey",linetype = "dashed",size=0.8) + 
  scale_x_log10() +
  theme_bw() 


uniCox_inhousedata <- uniCox
uniCox_inhousedata$dataset <- "inhouse_data"
save(uniCox_inhousedata, file = "d:/Work/Projects/BMC_Revised/Results/08.methods_compare/uniCox_inhousedata.RData")



# ============================================================================ #
# 3.roc ####
# ============================================================================ #
library(pROC)
df <- inhouse_merge

rocobj <- list()
cutOffPoint <- list()
cutOffPointText <- list()
auc_CI <- list()
for(method in colnames(inhouse_merge)[2:13]){
  # roc 分析
  rocobj[[method]] <- roc(df$Response, df[[method]], smooth = F)  # 曲线是否光滑，当光滑时，无法计算置信区间
  # 获取最大点用于后续绘图
  cutOffPoint[[method]] <- coords(rocobj[[method]], "best")
  cutOffPointText[[method]] <- paste0(round(cutOffPoint[[method]] [1],3),"(",round(cutOffPoint[[method]] [2],3),",",round(cutOffPoint[[method]] [3],3),")")
  # 计算 auc
  auc <- round(auc(rocobj[[method]])[1],3)
  auc_low <- round(ci(rocobj[[method]],of="auc")[1],3)
  auc_high <- round(ci(rocobj[[method]],of="auc")[3],3)
  auc_CI[[method]] <- paste0(auc," (",auc_low," - ",auc_high,")")
}

tmp <- as.data.frame(t(as.data.frame(auc_CI)))
tmp <- data.frame(method = rownames(tmp), roc = tmp$V1, dataset = "inhouse data")



ggroc(rocobj, aes = c("color"), size=1, legacy.axes = F)+ # FALSE时 横坐标为1-0 specificity；TRUE时 横坐标为0-1 1-specificity
  # 绘制对角线
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), colour='grey', linetype = 'dotdash') +
  # # 绘制临界点/阈值
  # geom_point(aes(x = cutOffPoint[[2]],y = cutOffPoint[[3]]))+ 
  # # geom_point(aes(x = cutOffPoint[[2]][[2]],y = cutOffPoint[[2]][[3]]))+
  # # geom_point(aes(x = cutOffPoint[[3]][[2]],y = cutOffPoint[[3]][[3]]))+
  # # geom_point(aes(x = cutOffPoint[[4]][[2]],y = cutOffPoint[[4]][[3]]))+
  # # 添加临界点/阈值文字标签
  # geom_text(aes(x = cutOffPoint[[2]],y = cutOffPoint[[3]],label=cutOffPointText),vjust=-1)+ 
  # # geom_text(aes(x = cutOffPoint[[2]][[2]],y = cutOffPoint[[2]][[3]],label=cutOffPointText[[2]]),vjust=-1)+ 
  # # geom_text(aes(x = cutOffPoint[[3]][[2]],y = cutOffPoint[[3]][[3]],label=cutOffPointText[[3]]),vjust=-1)+ 
  # # geom_text(aes(x = cutOffPoint[[4]][[2]],y = cutOffPoint[[4]][[3]],label=cutOffPointText[[4]]),vjust=-1)+ 
  # # # 添加 AUC
  # geom_text(aes(x = 0.3,y = 0.35,label = paste0("AUC of ",names(auc_CI), " = ", auc_CI)),vjust=-1)+ 
  # # geom_text(aes(x = 0.3,y = 0.25,label = paste0("AUC of ",names(auc_CI)[[2]], " = ", auc_CI[[2]])),vjust=-1)+ 
  # # geom_text(aes(x = 0.3,y = 0.15,label = paste0("AUC of ",names(auc_CI)[[3]], " = ", auc_CI[[3]])),vjust=-1)+ 
  # # geom_text(aes(x = 0.3,y = 0.05,label = paste0("AUC of ",names(auc_CI)[[4]], " = ", auc_CI[[4]])),vjust=-1)+ 
  # # 改配色
  # ggsci::scale_color_nejm()+
  theme_classic()

