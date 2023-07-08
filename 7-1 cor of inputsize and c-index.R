# input size and c-index

c_index_df <- openxlsx::read.xlsx("d:/Work/Projects/BMC_Revised/Results/07.tables/Table S3 HR and C-index of PScore and other pyroptosis-related models.xlsx")

unique(c_index_df$Method)
# [1] "PScore"                        "Wang_Ann.Transl.Med_2022"      "Shi_Medicine_2022"             "Niu_Math.Biosci.Eng_2022"     
# [5] "Ju_Front.Oncol_2021"           "Zhu_Stem.Cells.Int_2023"       "Xu_Front.Med_2022"             "Wu_PeerJ_2021"                
# [9] "Wu_Cancer.Med_2022"            "Li_Int.J.Gen.Med_2022"         "Wang_J.Investig.Dermatol_2022"

# PScore                     31
# Wang_Ann.Transl.Med_2022   6
# Shi_Medicine_2022:         3
# Niu_Math.Biosci.Eng_2022:  4
# Ju_Front.Oncol_2021:       9
# Zhu_Stem.Cells.Int_2023:   8
# Xu_Front.Med_2022:         8
# Wu_PeerJ_2021:             4
# Wu_Cancer.Med_2022:        4
# Li_Int.J.Gen.Med_2022:     5
# Wang_J.Investig.Dermatol_2022:  2


library(dplyr)
input_size <- data.frame("input_size" = c(31,6,3,4,9,8,8,4,4,5,2))
input_size <- apply(input_size, 1, function(.x)rep(.x,4)) %>% reshape::melt()
c_index_df$input_size <- input_size$value


library(ggplot2)
pdf("d:/Work/Projects/BMC_Revised/Results/06.figures/cor_inputsize_cindex.pdf",width = 4,height = 4)
ggplot(c_index_df, aes(x=input_size, y=Concordance))+
  geom_point()+
  ggpubr::stat_cor(data=c_index_df, method = "pearson")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()

