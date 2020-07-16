library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(ggpubr)

prot.abund <- read.table("../Data/Expression/PUF_time_course_hpt_LLCtherm_min_2_prot_log2_Imputed_remove_crap.tsv",sep="\t",header=T,row.names=1,stringsAsFactors = F)
pufs <- read.table("../Data/all_pufs_dufs.txt",sep="",header=F,stringsAsFactors = F)

t.prot.abund <- as.data.frame(t(prot.abund))
cols <- rownames(t.prot.abund)
sample.names <- unlist(strsplit(cols,split = "_rep[0-9]+"))
t.prot.abund["Sample.names"] <- sample.names

x <- t.prot.abund %>% group_by(Sample.names) %>% summarise_all(mean)
x <- as.data.frame(t(x))
colnames(x) <- as.character(unlist(unname(x[1,])))
x <- x[!(rownames(x) %in% c("Sample.names")),]
x[] <- lapply(x, function(i) as.numeric(as.character(i)))


x[rownames(x) %in% pufs[,1],"Protein"] <- "PUF"
x[!(rownames(x) %in% pufs[,1]),"Protein"] <- "Annotated"

x[,"Protein.Name"] <- rownames(x)

x_long <-x %>% gather(key="Sample",value="Abundance",hpt_T1,hpt_T2,hpt_T3,LL1210_T1,LL1210_T2,LL1210_T3,LL1210_T4)



p <- ggplot(x_long,aes(Sample,Abundance,fill=Protein))
p <- p + geom_boxplot(notch = T) + theme_cowplot() + ylim(c(17,37)) + stat_compare_means(aes(group=Protein),label="p.signif",label.y=35)
p <- p + ggtitle("Comparision of Protein Abundances:\nAnnotated vs. PUFs") + theme(axis.text.x=element_text(angle=45,vjust=0.5))
ggsave2(p,filename = "Images/annotated_vs_pufs_by_sample.png",dpi=300)



