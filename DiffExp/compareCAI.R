library(ggplot2)
library(cowplot)
library(ggpubr)

pufs <- read.table("../Data/all_pufs_dufs.txt",sep="",header=F,stringsAsFactors = F)
cai <- read.table("../Data/ctherm_cai.tsv",sep="\t",header=T,stringsAsFactors = F)

cai[cai[,1] %in% pufs[,1],"Type"] <- "PUF"
cai[!(cai[,1] %in% pufs[,1]),"Type"] <- "Annotated"

p <- ggplot(cai,aes(x=Type,y=cai))
p <- p + geom_boxplot(notch = T) + theme_cowplot() + stat_compare_means(label.y=0.80,label.x=1.5)
p <- p + ggtitle("Comparision of CAI:\nAnnotated vs. PUFs")
ggsave2(p,filename = "Images/annotated_vs_pufs_cai.png",dpi=300)


prot.abund <- read.table("../Data/Expression/PUF_time_course_hpt_LLCtherm_min_2_prot_log2_Imputed_remove_crap.tsv",sep="\t",header=T,row.names=1,stringsAsFactors = F)

t.prot.abund <- as.data.frame(t(prot.abund))
cols <- rownames(t.prot.abund)
sample.names <- unlist(strsplit(cols,split = "_rep[0-9]+"))
t.prot.abund["Sample.names"] <- sample.names

x <- t.prot.abund %>% group_by(Sample.names) %>% summarise_all(mean)
x <- as.data.frame(t(x))
colnames(x) <- as.character(unlist(unname(x[1,])))
x <- x[!(rownames(x) %in% c("Sample.names")),]
x[] <- lapply(x, function(i) as.numeric(as.character(i)))
x[,"Protein"] <- rownames(prot.abund)


cai <- cai[!(duplicated(cai[,1])),]

prot.cai <- merge(x,cai,by="Protein")

prot.cai.annot <- prot.cai[which(prot.cai$Type=="Annotated"),]
prot.cai.puf <- prot.cai[which(prot.cai$Type=="PUF"),]


pearson.corr.annot <- rep(0,6)
for (i in 2:7)
{
  tmp <- cor.test(prot.cai.annot[,i],prot.cai.annot$cai)
  pearson.corr.annot[i-1] <- tmp$estimate
}

pearson.corr.puf <- rep(0,6)
for (i in 2:7)
{
  tmp <- cor.test(prot.cai.puf[,i],prot.cai.puf$cai)
  pearson.corr.puf[i-1] <- tmp$estimate
}

corr.df <- data.frame(Correlation=c(pearson.corr.annot,pearson.corr.puf),Type=c(rep("Annotated",6),rep("PUF",6)),stringsAsFactors = F)

png("Images/ctherm_cai_vs_abundance_correlation_annotated_vs_puf.png")
boxplot(Correlation~Type,data=corr.df,col="blue",main="Comparing Correlation between CAI and Protein Abundance",ylab="Pearson Correlation",xlab="")
dev.off()


