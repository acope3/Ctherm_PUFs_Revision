library(ggfortify)

prot <- read.table("Expression/PUF_time_course_hpt_LLCtherm_min_2_prot_log2_Imputed_remove_crap.csv",sep=",",header=T,stringsAsFactors = F,row.names=1)
prot <- as.data.frame(t(prot),stringsAsFactors = F)
pca_res <- prcomp(prot,scale.=T)
prot[,"Sample"] <- c(rep("hpt_T1",4),rep("hpt_T2",4),rep("hpt_T3",4),rep("LL1210_T1",4),rep("LL1210_T2",4),rep("LL1210_T3",4),rep("LL1210_T4",4))
p<-autoplot(pca_res,data=prot,colour="Sample") + theme_bw()

ggsave(filename = "pca_abundancs.png",plot = p,dpi = 600)
