# Aim: Adjust the name of the gene set to match the TCGA-SKCM dataset


rm(list=ls())
path <- "d:/Work/Projects/BMC_Revised/"
mela_exper <- read.table("d:/Projects/Melanoma_Pyroptosis/Genelist/melanoma_mRNA.txt")  #FPKM
interest_genes <- openxlsx::read.xlsx(paste0(path, "Data/Gene_list/genelist_after_revise.xlsx"))  # genelist_after_revise.xlsx is Table 1

interest_genes <- interest_genes[1:75,]
interest_genes$Symbol.ID[interest_genes$Gene.ID %in% mela_exper$entrez_id == FALSE] 
# [1] "MIR223"
PRGs_74 <- interest_genes[interest_genes$Symbol.ID != "MIR223",]
# > dim(PRGs_74)
# [1] 74  5

save(PRGs_74, file = paste0(path, "Data/Gene_list/PRGs_74.RData"))
     