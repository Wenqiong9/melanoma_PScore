rm(list=ls())
load("d:/Work/Projects/BMC_Revised/Results/08.methods_compare/uniCox_GSE65904.RData")
load("d:/Work/Projects/BMC_Revised/Results/08.methods_compare/uniCox_GSE54467.RData")
load("d:/Work/Projects/BMC_Revised/Results/08.methods_compare/uniCox_TCGA.RData")
load("d:/Work/Projects/BMC_Revised/Results/08.methods_compare/uniCox_inhousedata.RData")

unicox_list <- list(GSE65904 = uniCox_GSE65904, inhousedata = uniCox_inhousedata,
                    GSE54467 = uniCox_GSE54467, TCGA = uniCox_TCGA)

uniCox <- data.frame()
for(i in names(unicox_list)){
  tmp <- unicox_list[[i]]
  tmp$gene[tmp$gene == "PRGs_prot"] <- "PScore" 
  tmp$method_dataset <- paste0(tmp$gene, "_", tmp$dataset)
  uniCox <- rbind(uniCox,tmp) 
}

uniCox$Risk_group <- factor(uniCox$Risk_group)
unicox_count <- data.frame()
for(i in unique(uniCox$gene)){
  tmp <- uniCox[uniCox$gene == i,]
  tmp$count <- table(tmp$Risk_group)[["Not_sig"]]
  unicox_count <- rbind(unicox_count, tmp)
}

order_methods <- rev(unique(unicox_count$gene[order(unicox_count$count)]))
unicox_order <- data.frame()
for(i in order_methods){
  tmp <- uniCox[uniCox$gene == i,]
  tmp <- tmp[order(tmp$Concordance),] 
  unicox_order <- rbind(unicox_order, tmp)
}
unicox_order$gene <- factor(unicox_order$gene, levels = order_methods)
unicox_order$method_dataset <- factor(unicox_order$method_dataset, levels = unicox_order$method_dataset)


library(ggplot2)
ggplot(unicox_order, aes(HR, method_dataset)) + 
  # facet_grid(~dataset) +
  geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
  geom_point(shape=18, aes(col=Risk_group, size = P_value))  +
  geom_vline(aes(xintercept = 1),color="darkgrey",linetype = "dashed",size=0.8) + 
  # scale_x_log10() +
  scale_color_manual(values = c("#FFDC91FF","#0072B5FF","#BC3C29FF"))+
  theme_bw() 


pdf("d:/Work/Projects/BMC_Revised/Results/08.methods_compare/forestplot2.pdf", width = 7, height = 8)
ggplot(unicox_order[1:22,], aes(HR, method_dataset)) +   ##  The x-axis label should also be log10 converted
  # facet_grid(~dataset) +
  geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
  geom_point(shape=18, aes(col=Risk_group, size = P_value))  +
  geom_vline(aes(xintercept = 1),color="darkgrey",linetype = "dashed",size=0.8) + 
  scale_x_log10() +
  labs(x = "log10(HR)") +
  scale_color_manual(values = c("#FFDC91FF","#0072B5FF","#BC3C29FF"))+
  theme_bw() 
ggplot(unicox_order[23:44,], aes(HR, method_dataset)) + 
  # facet_grid(~dataset) +
  geom_errorbarh(aes(xmax = HR.95L, xmin = HR.95H,col= Risk_group),height=0.3,size=0.7) +
  geom_point(shape=18, aes(col=Risk_group, size = P_value))  +
  geom_vline(aes(xintercept = 1),color="darkgrey",linetype = "dashed",size=0.8) + 
  scale_x_log10() +
  labs(x = "log10(HR)") +
  scale_color_manual(values = c("#FFDC91FF","#0072B5FF","#BC3C29FF"))+
  theme_bw() 
dev.off()


unicox_order$gene <- factor(unicox_order$gene, levels = rev(order_methods))
pdf("d:/Work/Projects/BMC_Revised/Results/08.methods_compare/barplot.pdf", width = 16, height = 6)
ggplot(unicox_order, aes(x=gene, y=Concordance, fill = dataset)) +
  geom_bar(stat="identity", position="dodge", width = .7)+
  geom_hline(yintercept=min(unicox_order$Concordance[unicox_order$gene == "PScore"])) +
  ggsci::scale_fill_nejm() +
  labs(x="")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, hjust=1))
dev.off()


colnames(unicox_order)[1] <- "Method"
unicox_order <- unicox_order[,-10]
unicox_order <- unicox_order[44:1,]
openxlsx::write.xlsx(unicox_order, "d:/Work/Projects/BMC_Revised/Results/07.tables/Table S3 HR and C-index of PScore and other pyroptosis-related models.xlsx")
