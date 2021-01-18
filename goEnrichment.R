library(topGO)
rm(list = ls())

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

fisherTestForClust <- function()
{
  ## Change this to BP for biological processes instead of molecular function
  ontology <- "BP"
  
  pufs <- read.table("Data/measured_pufs.txt",sep="",header=F,stringsAsFactors = F)
  pufs <- pufs[,1]
  
  geneID2GO <- readMappings(file="Data/GO/ctherm_protein_to_go_map_updated.tsv")
  
  prot <- read.table("WGCNA/Cleaned_data/hpt_cleaned_imputed_data_transposed.tsv",sep="\t",header=T,stringsAsFactors = F,row.names=1)
  gene.universe <- rownames(prot)
  
  clusters <- read.table("Clust/hpt_no_norm/Clusters_Objects.tsv",sep="\t",header = T,skip = 1,stringsAsFactors = F)
  
  for (i in 1:ncol(clusters))
  {
    interesting.genes <- clusters[,i]
    
    x <- factor(as.integer(gene.universe %in% interesting.genes))
    names(x) <- gene.universe
    
    GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=10)
    elimResultFisher <- runTest(GOdata,statistic = "fisher",algorithm = "elim")
    weight01ResultFisher<- runTest(GOdata,statistic = "fisher",algorithm = "weight01")
    allRes <- GenTable(GOdata,weight01Fisher=weight01ResultFisher,elimFisher=elimResultFisher,topNodes=100,orderBy="weight01Fisher")
    sig.go.terms <- allRes[which(allRes$weight01Fisher < 0.05),1]
    if (length(sig.go.terms) > 0)
    {
      for (j in 1:length(sig.go.terms))
      {
        genes <- genesInTerm(GOdata,sig.go.terms[j])[[1]]
        genes.int <- genes[which(genes %in% interesting.genes)]
        genes.tog <- paste0(genes.int,collapse=";")
        allRes[j,"Protein"] <- genes.tog
        y<-prot[unlist(genes.int),]
        pufs.str <- paste0(pufs[which(pufs %in% interesting.genes)],collapse=";")
        allRes[j,"PUF"] <- pufs.str
      }
    }
    write.table(allRes,file.path("Clust","hpt_no_norm","GO",paste0("go_enrichment_",ontology,"_cluster_",i-1,".tsv")),sep="\t",col.names = T,row.names=F,quote=F)
  }
}



topDiffGenes <- function(allScores)
{
  return(allScores < 0.05)
}

ksTestForDiffExp <- function(ontology,fold.change.file,output,gene.universe)
{
  geneID2GO <- readMappings(file="Data/GO/ctherm_protein_to_go_map_updated.tsv")
  
  fold.change <- read.table(fold.change.file,sep="\t",header = T,row.names=1,stringsAsFactors = F)
  gene.universe <- rownames(fold.change)
  
  x <- fold.change$adj.P.Val
  names(x) <- gene.universe
  
  GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=5,geneSelectionFun=topDiffGenes)
  classicResultKS <- runTest(GOdata,statistic = "ks",algorithm = "classic")
  elimResultKS <- runTest(GOdata,statistic = "ks",algorithm = "elim")
  weight01ResultKS <- runTest(GOdata,statistic = "ks",algorithm = "weight01")
  allRes <- GenTable(GOdata,weight01KS=weight01ResultKS,topNodes=100)
  sig.go.terms <- allRes[which(allRes$weight01KS < 0.05),1]
  for (i in 1:length(sig.go.terms))
  {
    genes <- genesInTerm(GOdata,sig.go.terms[i])[[1]]
    genes.tog <- paste0(genes,collapse=";")
    allRes[i,"Protein"] <- genes.tog
    x<-fold.change[unlist(genes),]
    mean.fc.all <- mean(x$logFC)
    mean.fc.puf <- mean(x[which(x$PUF == "Yes"),"logFC"])
    pufs <- paste0(rownames(x)[which(x$PUF == "Yes")],collapse=";")
    allRes[i,"PUF"] <- pufs
    allRes[i,"Mean.logFC.all"] <- mean.fc.all
    allRes[i,"Mean.logFC.puf"] <- mean.fc.puf
  }

  write.table(allRes,output,sep="\t",col.names = T,row.names=F,quote=F)
}

ksTestForDiffCoExp <- function()
{
  ontology <- "MF"
  
  geneID2GO <- readMappings(file="Data/GO/ctherm_protein_to_go_map_updated.tsv")
  
  fold.change <- read.table("DiffCoExpress/differential_coexpression_strains.tsv",sep="\t",row.names = 1,header = T,stringsAsFactors = F)
  fold.change <- fold.change[which(!is.na(fold.change$dC)),]
  gene.universe <- rownames(fold.change)

  x <- fold.change$q.value
  names(x) <- gene.universe
  
  GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=5,geneSelectionFun=topDiffGenes)
  classicResultFisher <- runTest(GOdata,statistic="fisher",algorithm="classic")
  elimResultFisher <- runTest(GOdata,statistic = "fisher",algorithm = "elim")
  classicResultKS <- runTest(GOdata,statistic = "ks",algorithm = "classic")
  elimResultKS <- runTest(GOdata,statistic = "ks",algorithm = "elim")
  weight01ResultKS <- runTest(GOdata,statistic = "ks",algorithm = "weight01")
  allRes <- GenTable(GOdata,classicFisher=classicResultFisher,elimFisher=elimResultFisher,classicKS=classicResultKS,elimKS=elimResultKS,weight01KS=weight01ResultKS,topNodes=100,orderBy="weight01KS")
  sig.go.terms <- allRes[which(allRes$weight01KS < 0.05),1]
  for (i in 1:length(sig.go.terms))
  {
    genes <- genesInTerm(GOdata,sig.go.terms[i])[[1]]
    genes.tog <- paste0(genes,collapse=";")
    allRes[i,"Protein"] <- genes.tog
    x<-fold.change[unlist(genes),]
    mean.fc.all <- mean(x$dC)
    mean.fc.puf <- mean(x[which(x$PUF == "Yes"),"dC"])
    pufs <- paste0(rownames(x)[which(x$PUF == "Yes")],collapse=";")
    allRes[i,"PUF"] <- pufs
    allRes[i,"Mean.dC.all"] <- mean.fc.all
    allRes[i,"Mean.dC.puf"] <- mean.fc.puf
  }
  
  write.table(allRes,file.path("DiffCoExpress/",paste0("go_enrichment_",ontology,".tsv")),sep="\t",col.names = T,row.names=F,quote=F)
}

ksTestForDiffExp("BP","DiffExp/Tables/hpt_vs_LL1210_early_log.tsv",file.path("DiffExp","GO","EL",paste0("go_enrichment_early_log_BP.tsv")))
ksTestForDiffExp("MF","DiffExp/Tables/hpt_vs_LL1210_early_log.tsv",file.path("DiffExp","GO","EL",paste0("go_enrichment_early_log_MF.tsv")))
ksTestForDiffExp("BP","DiffExp/Tables/hpt_vs_LL1210_mid_log.tsv",file.path("DiffExp","GO","ML",paste0("go_enrichment_mid_log_BP.tsv")))
ksTestForDiffExp("MF","DiffExp/Tables/hpt_vs_LL1210_mid_log.tsv",file.path("DiffExp","GO","ML",paste0("go_enrichment_mid_log_MF.tsv")))
ksTestForDiffExp("BP","DiffExp/Tables/hpt_vs_LL1210_late_log.tsv",file.path("DiffExp","GO","LL",paste0("go_enrichment_late_log_BP.tsv")))
ksTestForDiffExp("MF","DiffExp/Tables/hpt_vs_LL1210_late_log.tsv",file.path("DiffExp","GO","LL",paste0("go_enrichment_late_log_MF.tsv")))



#ksTestForDiffCoExp()
#fisherTestForClust()
#fisherTestForWGCNA()

# ontology <- "MF"
# 
# geneID2GO <- readMappings(file="Data/GO/ctherm_protein_to_go_map_updated.tsv")
# 
# interesting.genes <- read.table("Evolutionary/phylo_profiling_significant_correlated_proteins.txt",sep="",header = T,stringsAsFactors = F)
# interesting.genes <- interesting.genes[,1]
# all.prot <- read.table("all_protein_ids.txt",sep="",header=F,stringsAsFactors = F)
# gene.universe <- all.prot[,1]
# 
# x <- factor(as.integer(gene.universe %in% interesting.genes))
# names(x) <- gene.universe
# 
# 
# GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=5)
# classicResultFisher <- runTest(GOdata,statistic="fisher",algorithm="classic")
# elimResultFisher <- runTest(GOdata,statistic = "fisher",algorithm = "elim")
# weight01ResultFisher <- runTest(GOdata,statistic = "fisher",algorithm = "weight01")
# allRes <- GenTable(GOdata,classicFisher=classicResultFisher,elimFisher=elimResultFisher,weight01Fisher=weight01ResultFisher,topNodes=100,orderBy="weight01Fisher")
# write.table(allRes,file.path("Evolutionary/",paste0("go_enrichment_",ontology,".tsv")),sep="\t",col.names = T,row.names=F,quote=F)


