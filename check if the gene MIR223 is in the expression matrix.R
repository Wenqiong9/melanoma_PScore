# Aim: check if the gene MIR223 (ENSG00000284567) is in the expression matrix
# the expression matrix was downloaded from https://xenabrowser.net/datapages/ named 'TCGA Melanoma (SKCM) (17 datasets)'

rm(list=ls())
path <- "d:/Work/Projects/BMC_Revised_2/"
mela_exper <- read.table("d:/Work/Projects/BMC_Revised_2/TCGA-SKCM.htseq_fpkm.tsv.gz", header = TRUE)  #FPKM
head(mela_exper)[,1:4]

ensembl <- gsub("\\..$", "", mela_exper$Ensembl_ID)
ensembl <- gsub("\\...$", "", ensembl)

"ENSG00000284567" %in% ensembl
