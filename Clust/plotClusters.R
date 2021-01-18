library(ggplot2)
library(reshape2)
library(topGO)
library(cowplot)

getProtAssociatedWGo <- function(go.term,go.map)
{
  prot.has.go <- unlist(lapply(go.map,function(x){go.term %in% x}))
  return(names(go.map[prot.has.go]))
}

createGOdataObject<- function(geneList,id2Go,ontology="BP",...)
{
  GOData <- new("topGOdata",ontology=ontology,allGenes=geneList,gene2GO=id2Go,annot=annFUN.gene2GO,...)
  return(GOData)
}

createClusterPlot <- function(prot.abund,main.target=NULL,other.target=NULL,other.target.label=NULL)
{
  color.scheme <- c("grey","blue","red","green")
  if (is.null(main.target) && is.null(other.target))
  {
    color.scheme <- color.scheme[1:2]
    names(color.scheme) <- c("Annotated","PUF")
  } else if (is.null(other.target)) {
    color.scheme <- color.scheme[1:3]
    names(color.scheme) <- c("Annotated","PUF",main.target)
  } else{
    color.scheme <- color.scheme[1:4]
    names(color.scheme) <- c("Annotated","PUF",main.target,other.target.label)
  }
  
  if (!is.null(main.target)) 
  {
    prot.abund[main.target,"PUF"] <- main.target
    main.target.index <- which(prot.abund$Protein == main.target)
    puf.index <- which(prot.abund$PUF == "PUF" & prot.abund$Protein != main.target)
  } else{
    main.target.index <- c()
    puf.index <- which(prot.abund$PUF == "PUF")
  }
  
  if (!is.null(other.target))
  {
    prot.abund[other.target,"PUF"] <- other.target.label
    prot.abund[puf.index,"PUF"] <- "PUF"
    prot.abund[main.target.index,"PUF"] <- main.target
    go.index <- which(prot.abund$PUF == other.target.label)
  } else{
    go.index <- c()
  }
  prot.id <- prot.abund$Protein
  annot.index <- which(prot.abund$PUF == "Annotated")
  prot.abund[,"Protein"] <- factor(prot.abund$Protein,levels=c(prot.id[annot.index],prot.id[puf.index],prot.id[go.index],prot.id[main.target.index])) # make sure target puf is on top
  prot.abund.melt <- melt(prot.abund,id.vars=c("Protein","PUF"),variable.name="TimePoint",value.name="Abundance")
  p <- ggplot(prot.abund.melt,aes(x=TimePoint,y=Abundance,group=Protein,colour=PUF)) + geom_line() 
  p <- p + xlab("Time Point") + ylab("Normalized Abundance") + theme_bw()
  #p <- p + scale_color_manual(values=ifelse(is.null(main.target),c("Annotated"="grey","PUF"="blue"),c("Annotated"="grey","PUF"="blue",main.target = "red")))
  p <- p + scale_color_manual(values=color.scheme,name="Annotation")
  p <- p + theme(axis.text = element_text(size=8,face="bold",colour = "black"),axis.title = element_text(size=10,face="bold"),legend.title = element_text(size=12,face="bold"),legend.text = element_text(size=8,face="bold"))
  return(p)
  
}

produceSpecifiFiguresForPaper <- function()
{
  ontology.bp <- "BP"
  ontology.mf <- "MF"
  geneID2GO <- readMappings(file="../Data/GO/ctherm_protein_to_go_map_updated.tsv")
  puf <- read.table("../Data/measured_pufs.txt",sep="",header=F,stringsAsFactors = F)
  
  prot.hpt <- read.table("../WGCNA/Cleaned_data/hpt_cleaned_imputed_data_transposed.tsv",sep="\t",header=T,stringsAsFactors = F,row.names=1)
  gene.universe.hpt <- rownames(prot.hpt)
  prot.file <- "hpt/Processed_Data/hpt_cleaned_imputed_data_transposed.tsv_processed.tsv"
  clust.file <- "hpt/Clusters_Objects.tsv"
  prot.hpt.df <- read.table(prot.file,sep="\t",header=T,stringsAsFactors = F,row.names=1)
  clust.hpt.df <- read.table(clust.file,sep="\t",header=T,stringsAsFactors = F,skip = 1)
  prot.hpt.df[,"Protein"] <- row.names(prot.hpt.df)
  prot.hpt.df[,"PUF"] <- ifelse(prot.hpt.df$Protein %in% puf[,1],"PUF","Annotated")
  
  
  prot.ll <- read.table("../WGCNA/Cleaned_data/ll1210_cleaned_imputed_data_transposed.tsv",sep="\t",header=T,stringsAsFactors = F,row.names=1)
  gene.universe.ll <- rownames(prot.ll)
  prot.file <- "LL1210/Processed_Data/ll1210_cleaned_imputed_data_transposed.tsv_processed.tsv"
  clust.file <- "LL1210/Clusters_Objects.tsv"
  prot.ll.df <- read.table(prot.file,sep="\t",header=T,stringsAsFactors = F,row.names=1)
  clust.ll.df <- read.table(clust.file,sep="\t",header=T,stringsAsFactors = F,skip = 1)
  prot.ll.df[,"Protein"] <- row.names(prot.ll.df)
  prot.ll.df[,"PUF"] <- ifelse(prot.ll.df$Protein %in% puf[,1],"PUF","Annotated")
  
  
  ## WP_003512015.1
  interesting.genes.hpt <- clust.hpt.df[which(clust.hpt.df[,6] != ""),6]
  x.hpt <- factor(as.integer(gene.universe.hpt %in% interesting.genes.hpt))
  names(x.hpt) <- gene.universe.hpt
  
  GOdata.hpt <- createGOdataObject(x.hpt,geneID2GO,ontology = ontology.mf,nodeSize=10)
  other.genes <- genesInTerm(GOdata.hpt,"GO:0051539")[[1]]
  other.genes <- other.genes[which(other.genes %in% interesting.genes.hpt)]
  
  abund <- prot.hpt.df[interesting.genes.hpt,]
  p <- createClusterPlot(abund,main.target = "WP_003512015.1",other.target = other.genes,other.target.label="GO:0051539")  
  ggsave2("WP_003512015.1_hpt_clust.png",p,dpi=600)
  
  ## WP_003516357.1
  interesting.genes.ll <- clust.ll.df[which(clust.ll.df[,14] != ""),14]
  
  x.ll<- factor(as.integer(gene.universe.ll %in% interesting.genes.ll))
  names(x.ll) <- gene.universe.ll
  
  GOdata.ll <- createGOdataObject(x.ll,geneID2GO,ontology = ontology.bp,nodeSize=10)
  other.genes <- genesInTerm(GOdata.ll,"GO:0000272")[[1]]
  other.genes <- other.genes[which(other.genes %in% interesting.genes.ll)]
  abund <- prot.ll.df[interesting.genes.ll,]
  p <- createClusterPlot(abund,main.target = "WP_003516357.1",other.target = other.genes,other.target.label="GO:0000272")  
  ggsave2("WP_003516357.1_ll1210_clust.png",p,dpi=600)
  
  
  ## WP_003511984.1
  interesting.genes.hpt <- clust.hpt.df[which(clust.hpt.df[,1] != ""),1]
  
  x.hpt <- factor(as.integer(gene.universe.hpt %in% interesting.genes.hpt))
  names(x.hpt) <- gene.universe.hpt
  
  GOdata.hpt <- createGOdataObject(x.hpt,geneID2GO,ontology = ontology.bp,nodeSize=10)
  other.genes <- genesInTerm(GOdata.hpt,"GO:0009057")[[1]]
  other.genes <- other.genes[which(other.genes %in% interesting.genes.hpt)]
  abund <- prot.hpt.df[interesting.genes.hpt,]
  p <- createClusterPlot(abund,main.target = "WP_003511984.1",other.target = other.genes,other.target.label="GO:0009057")  
  ggsave2("WP_003511984.1_hpt_clust.png",p,dpi=600)
  
  
  ## WP_003519433.1
  interesting.genes.ll <- clust.ll.df[which(clust.ll.df[,2] != ""),2]
  
  x.ll<- factor(as.integer(gene.universe.ll %in% interesting.genes.ll))
  names(x.ll) <- gene.universe.ll
  
  GOdata.ll <- createGOdataObject(x.ll,geneID2GO,ontology = ontology.mf,nodeSize=10)
  other.genes <- genesInTerm(GOdata.ll,"GO:0016740")[[1]]
  other.genes <- other.genes[which(other.genes %in% interesting.genes.ll)]
  abund <- prot.ll.df[interesting.genes.ll,]
  p <- createClusterPlot(abund,main.target = "WP_003519433.1",other.target = other.genes,other.target.label="GO:0016740")  
  ggsave2("WP_003519433.1_ll1210_clust.png",p,dpi=600)
  
}

produceSpecifiFiguresForPaper()

# prot <- read.table("../WGCNA/Cleaned_data/hpt_cleaned_imputed_data_transposed.tsv",sep="\t",header=T,stringsAsFactors = F,row.names=1)
# gene.universe <- rownames(prot)
# ontology <- "MF"
# geneID2GO <- readMappings(file="../Data/GO/ctherm_protein_to_go_map_updated.tsv")
# 
# prot.file <- "hpt/Processed_Data/hpt_cleaned_imputed_data_transposed.tsv_processed.tsv"
# clust.file <- "hpt/Clusters_Objects.tsv"
# prot.df <- read.table(prot.file,sep="\t",header=T,stringsAsFactors = F,row.names=1)
# clust.df <- read.table(clust.file,sep="\t",header=T,stringsAsFactors = F,skip = 1)
# 
# 
# prot.df[,"Protein"] <- row.names(prot.df)
# prot.df[,"PUF"] <- ifelse(prot.df$Protein %in% puf[,1],"PUF","Annotated")
# 
# 
# interesting.genes <- clust.df[which(clust.df[,6] != ""),6]
# 
# 
# 
# x <- factor(as.integer(gene.universe %in% interesting.genes))
# names(x) <- gene.universe
# 
# GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=10)
# other.genes <- genesInTerm(GOdata,"GO:0051539")[[1]]
# other.genes <- other.genes[which(other.genes %in% interesting.genes)]
# 
# clust.plot <- vector(mode="list",length=ncol(clust.df))
# for (i in 1:ncol(clust.df))
# {
#   abund <- prot.df[clust.df[which(clust.df[,i] != ""),i],]
#   clust.plot[[i]] <- createClusterPlot(abund)  
# }
# 
# legend <- get_legend(clust.plot[[1]]+ guides(colour = guide_legend(override.aes = list(size=10))) + theme(legend.title=element_text(size=20),legend.text=element_text(size=14)))
# 
# top <- plot_grid(clust.plot[[1]] + ggtitle("hpt_C0") + theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()) ,
#                  clust.plot[[2]] + ggtitle("hpt_C1") + theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()) ,
#                  clust.plot[[3]] + ggtitle("hpt_C2") + theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()) ,
#                  clust.plot[[4]] + ggtitle("hpt_C3")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()),
#                  labels=c("A","B","C","D"),nrow=1)
# mid <- plot_grid(clust.plot[[5]]+ ggtitle("hpt_C4") + theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()),
#                  clust.plot[[6]]+ ggtitle("hpt_C5")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()),
#                  clust.plot[[7]]+ ggtitle("hpt_C6")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()),
#                  clust.plot[[8]]+ ggtitle("hpt_C7")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()),
#                  labels=c("E","F","G","H"),nrow = 1)
# bottom <- plot_grid(clust.plot[[9]]+ ggtitle("hpt_C8")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()),
#                     clust.plot[[10]]+ ggtitle("hpt_C9")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()),
#                     clust.plot[[11]]+ ggtitle("hpt_C10")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.title.x=element_blank()),
#                     legend,
#                     labels=c("I","J","K",""),nrow=1)
# 
# all <- plot_grid(top,mid,bottom,ncol = 1)
# 
# ggsave2("hpt_clusters.png",plot = all,dpi=600,height=12,width=12)
# 
# 
# prot.file <- "LL1210/Processed_Data/ll1210_cleaned_imputed_data_transposed.tsv_processed.tsv"
# clust.file <- "LL1210/Clusters_Objects.tsv"
# prot.df <- read.table(prot.file,sep="\t",header=T,stringsAsFactors = F,row.names=1)
# clust.df <- read.table(clust.file,sep="\t",header=T,stringsAsFactors = F,skip = 1)
# prot.df[,"Protein"] <- row.names(prot.df)
# prot.df[,"PUF"] <- ifelse(prot.df$Protein %in% puf[,1],"PUF","Annotated")
# 
# clust.plot <- vector(mode="list",length=ncol(clust.df))
# for (i in 1:ncol(clust.df))
# {
#   abund <- prot.df[clust.df[which(clust.df[,i] != ""),i],]
#   clust.plot[[i]] <- createClusterPlot(abund)  
# }
# 
# legend <- get_legend(clust.plot[[1]]+ guides(colour = guide_legend(override.aes = list(size=10))) + theme(legend.title=element_text(size=20),legend.text=element_text(size=14)))
# 
# top <- plot_grid(clust.plot[[1]] + ggtitle("LL1210_C0") + theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()) ,
#                  clust.plot[[2]] + ggtitle("LL1210_C1") + theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()) ,
#                  clust.plot[[3]] + ggtitle("LL1210_C2") + theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()) ,
#                  clust.plot[[4]] + ggtitle("LL1210_C3")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                  labels=c("A","B","C","D"),nrow=1)
# mid.1 <- plot_grid(clust.plot[[5]]+ ggtitle("LL1210_C4") + theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                  clust.plot[[6]]+ ggtitle("LL1210_C5")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                  clust.plot[[7]]+ ggtitle("LL1210_C6")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                  clust.plot[[8]]+ ggtitle("LL1210_C7")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                  labels=c("E","F","G","H"),nrow=1)
# mid.2 <- plot_grid(clust.plot[[9]]+ ggtitle("LL1210_C8") + theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                    clust.plot[[10]]+ ggtitle("LL1210_C9")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                    clust.plot[[11]]+ ggtitle("LL1210_C10")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                    clust.plot[[12]]+ ggtitle("LL1210_C11")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                    labels=c("I","J","K","L"),nrow=1)
# bottom <- plot_grid(clust.plot[[13]]+ ggtitle("LL1210_C12")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                     clust.plot[[14]]+ ggtitle("LL1210_C13")+ theme(legend.position = "none",axis.title=element_text(size=6),axis.text=element_text(size=8),axis.title.x=element_blank()),
#                     legend,
#                     labels=c("M","N",""),nrow=1)
# 
# all <- plot_grid(top,mid.1,mid.2,bottom,ncol = 1)
# 
# ggsave2("ll1210_clusters.png",plot = all,dpi=600,height=12,width=12)



