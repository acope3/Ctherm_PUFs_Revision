library(topGO)
rm(list = ls())

createGOdataObject<- function(geneList,id2Go,ontology="BP",...)
{
  GOData <- new("topGOdata",ontology=ontology,allGenes=geneList,gene2GO=id2Go,annot=annFUN.gene2GO,...)
  return(GOData)
}

fisherTestForClust <- function()
{
  ## Change this to BP for biological processes instead of molecular function
  ontology <- "MF"
  
  geneID2GO <- readMappings(file="Data/GO/ctherm_protein_to_go_map_updated.tsv")
  
  prot <- read.table("WGCNA/Cleaned_data/hpt_cleaned_imputed_data_transposed.tsv",sep="\t",header=T,stringsAsFactors = F)
  gene.universe <- prot$Protein
  
  clusters <- read.table("Clust/hpt/Clusters_Objects.tsv",sep="\t",header = T,skip = 1)
  
  for (i in 1:ncol(clusters))
  {
    interesting.genes <- clusters[,i]
    
    x <- factor(as.integer(gene.universe %in% interesting.genes))
    names(x) <- gene.universe
    
    GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=10)
    resultFisher <- runTest(GOdata,statistic = "fisher",algorithm = "elim")
    allRes <- GenTable(GOdata,Fisher=resultFisher,topNodes=100)
    write.table(allRes,file.path("Clust","hpt","GO",paste0("go_enrichment_",ontology,"_cluster_",i-1,".tsv")),sep="\t",col.names = T,row.names=F,quote=F)
  }
}

fisherTestForWGCNA <- function()
{
  ontology <- "BP"
  
  geneID2GO <- readMappings(file="Data/GO/ctherm_protein_to_go_map_updated.tsv")
  
  prot <- read.table("WGCNA/Cleaned_data/ll1210_cleaned_imputed_data_transposed.tsv",sep="\t",header=T,stringsAsFactors = F)
  gene.universe <- prot$Protein
  
  clusters <- read.table("WGCNA/Modules/LL1210/ll1210_modules.tsv",sep="\t",header = T)
  
  for (i in unique(clusters$Cluster))
  {
    interesting.genes <- clusters[which(clusters$Cluster==i),"Protein"]
    print(interesting.genes)
    x <- factor(as.integer(gene.universe %in% interesting.genes))
    names(x) <- gene.universe
    
    GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=10)
    resultFisher <- runTest(GOdata,statistic = "fisher",algorithm = "elim")
    allRes <- GenTable(GOdata,Fisher=resultFisher,topNodes=100)
    write.table(allRes,file.path("WGCNA","Modules","LL1210","GO",paste0("go_enrichment_",ontology,"_cluster_",i,".tsv")),sep="\t",col.names = T,row.names=F,quote=F)
  }
}


topDiffGenes <- function(allScores)
{
  return(allScores < 0.05)
}

ksTestForDiffExp <- function()
{
  geneID2GO <- readMappings(file="Data/GO/ctherm_protein_to_go_map_updated.tsv")
  
  fold.change <- read.table("DiffExp/Tables/hpt_vs_LL1210_late_log.tsv",sep="\t",header = T,stringsAsFactors = F)
  gene.universe <- fold.change$Protein
  
  x <- fold.change$adj.P.Val
  names(x) <- gene.universe
  
  GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=5,geneSelectionFun=topDiffGenes)
  classicResultKS <- runTest(GOdata,statistic = "ks",algorithm = "classic")
  elimResultKS <- runTest(GOdata,statistic = "ks",algorithm = "elim")
  weight01ResultKS <- runTest(GOdata,statistic = "ks",algorithm = "weight01")
  allRes <- GenTable(GOdata,classicKS=classicResultKS,elimKS=elimResultKS,weight01KS=weight01ResultKS,topNodes=100)
  write.table(allRes,file.path("DiffExp",paste0("go_enrichment_late_log_",ontology,".tsv")),sep="\t",col.names = T,row.names=F,quote=F)
}

ksTestForDiffCoExp <- function()
{
  ontology <- "MF"
  
  geneID2GO <- readMappings(file="Data/GO/ctherm_protein_to_go_map_updated.tsv")
  
  fold.change <- read.table("DiffCoExpress/differential_coexpression_strains.tsv",sep="\t",header = T,stringsAsFactors = F)
  gene.universe <- fold.change$Protein
  
  
  x <- fold.change$q.value
  names(x) <- gene.universe
  
  GOdata <- createGOdataObject(x,geneID2GO,ontology = ontology,nodeSize=5,geneSelectionFun=topDiffGenes)
  classicResultFisher <- runTest(GOdata,statistic="fisher",algorithm="classic")
  elimResultFisher <- runTest(GOdata,statistic = "fisher",algorithm = "elim")
  classicResultKS <- runTest(GOdata,statistic = "ks",algorithm = "classic")
  elimResultKS <- runTest(GOdata,statistic = "ks",algorithm = "elim")
  weight01ResultKS <- runTest(GOdata,statistic = "ks",algorithm = "weight01")
  allRes <- GenTable(GOdata,classicFisher=classicResultFisher,elimFisher=elimResultFisher,classicKS=classicResultKS,elimKS=elimResultKS,weight01KS=weight01ResultKS,topNodes=100,orderBy="weight01KS")
  write.table(allRes,file.path("DiffCoExpress/",paste0("go_enrichment_",ontology,".tsv")),sep="\t",col.names = T,row.names=F,quote=F)
}

fisherTestForClust()
