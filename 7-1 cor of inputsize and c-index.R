# input size and c-index

c_index_df <- openxlsx::read.xlsx("d:/Work/Projects/BMC_Revised/Results/07.tables/Table S3 HR and C-index of PScore and other pyroptosis-related models.xlsx")

unique(c_index_df$Method)

# Wang_J.Investig.Dermatol_2022:  2
# Li_Int.J.Gen.Med_2022:     5
# Wu_Cancer.Med_2022:        4
# Ju_Front.Oncol_2021:       9
# Wu_PeerJ_2021:             4
# Niu_Math.Biosci.Eng_2022:  4
# Xu_Front.Med_2022:         8
# Zhu_Stem.Cells.Int_2023:   8
# Shi_Medicine_2022:         3
# Wang_Ann.Transl.Med_2022   6
# PScore                     31
# PRGs_74                    74

iput_size <- data.frame("iput_size" = c(2,5,4,9,4,4,8,8,3,6,31,74))
iput_size <- apply(iput_size, 1, function(.x)rep(.x,4)) %>% reshape::melt()
c_index_df$iput_size <- iput_size$value


library(ggplot2)
ggplot(c_index_df, aes(x=iput_size, y=Concordance))+
  geom_point()+
  # xlab("Pyroptosis Score")+
  # ylab(i)+
  stat_smooth(method = "lm",se=T)+
  ggpubr::stat_cor(data=c_index_df, method = "pearson")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())


