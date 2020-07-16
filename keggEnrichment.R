library(clusterProfiler)

kegg_overRepresentation <- function(data,module=F)
{
  if (module)
  {
    enrichment <- enrichMKEGG(data,organism = "ctx",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1)
  } else{
    enrichment <- enrichKEGG(data,organism = "ctx",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.1)
  }
  return(enrichment)
}

diffExp.gsea <- function(diff.exp.file,target.directory,module=F)
{
  data <- read.table(diff.exp.file,sep="\t",header=T,stringsAsFactors = F)
  old_ids <- read.table("old_id_map.tsv",sep="\t",header=T)
  new.data <- merge(data,old_ids,how="left",on="Protein")
  new.data <- new.data[-which(duplicated(new.data$Locus)),]
  fold.change <- new.data$logFC
  names(fold.change) <- new.data$Locus
  fold.change <- sort(fold.change,decreasing = T)
  if (module)
  {
    gsea <- gseMKEGG(fold.change,organism="ctx")
  } else{
    gsea <- gseKEGG(fold.change,organism = "ctx")
  }
  write.table(gsea@result,file.path(target.directory,ifelse(module,"kegg_gsea_module.tsv","kegg_gsea.tsv")),sep = "\t",col.names = T,row.names = F,quote=F)
}

diffExp.enrichment <- function(diff.exp.file,target.directory,module=F)
{
  data <- read.table(diff.exp.file,sep="\t",header=T,stringsAsFactors = F)
  old_ids <- read.table("old_id_map.tsv",sep="\t",header=T)
  new.data <- merge(data,old_ids,how="left",on="Protein")
  diff.exp <- new.data[which(new.data$adj.P.Val < 0.05),"Locus"]
  enriched.kegg <- kegg_overRepresentation(diff.exp,module)
  if(!is.null(enriched.kegg))
  {
    write.table(enriched.kegg@result,file.path(target.directory,paste0(ifelse(module,"kegg_enrichment_module.tsv","kegg_enrichment.tsv"))),sep="\t",row.names = F,col.names = T,quote=F)
  }
}

coExpress.clust.enrichment <- function(clust.file,target.directory,module=F)
{
  data <- read.table(clust.file,sep="\t",header=T,stringsAsFactors = F,skip=1)
  old_ids <- read.table("old_id_map.tsv",sep="\t",header=T)
  for (i in 1:ncol(data))
  {
    new.data <- data[,i]
    interesting.genes <- old_ids[which(old_ids$Protein %in% new.data),"Locus"]
    enriched.kegg <- kegg_overRepresentation(interesting.genes,module)
    if(!is.null(enriched.kegg))
    {
      write.table(enriched.kegg@result,file.path(target.directory,paste0(ifelse(module,"kegg_enrichment_module_clust_","kegg_enrichment_clust_"),i-1,".tsv")),sep="\t",row.names = F,col.names = T,quote=F)
    }
  }
}

coExpress.WGCNA.enrichment <- function(clust.file,target.directory,module=F)
{
  data <- read.table(clust.file,sep="\t",header=T,stringsAsFactors = F)
  old_ids <- read.table("old_id_map.tsv",sep="\t",header=T)
  for (i in unique(data$Cluster))
  {
    new.data <- data[which(data$Cluster == i),"Protein"]
    interesting.genes <- old_ids[which(old_ids$Protein %in% new.data),"Locus"]
    enriched.kegg <- kegg_overRepresentation(interesting.genes,module)
    if(!is.null(enriched.kegg))
    {
      write.table(enriched.kegg@result,file.path(target.directory,paste0(ifelse(module,"kegg_enrichment_module_WGCNA_","kegg_enrichment_WGCNA_"),i,".tsv")),sep="\t",row.names = F,col.names = T,quote=F)
    }
  }
}


diffCoExpress.enrichment <- function(diff.exp.file,target.directory,module=F)
{
  data <- read.table(diff.exp.file,sep="\t",header=T,stringsAsFactors = F)
  old_ids <- read.table("old_id_map.tsv",sep="\t",header=T)
  new.data <- merge(data,old_ids,how="left",on="Protein")
  diff.exp <- new.data[which(new.data$q.value < 0.05),"Locus"]
  enriched.kegg <- kegg_overRepresentation(diff.exp,module)
  write.table(enriched.kegg@result,file.path(target.directory,ifelse(module,"kegg_enrichment_module.tsv","kegg_enrichment.tsv")),sep="\t",row.names = F,col.names = T,quote=F)
}

diffCoExpress.gsea <- function(diff.exp.file,target.directory,module=F)
{
  data <- read.table(diff.exp.file,sep="\t",header=T,stringsAsFactors = F)
  old_ids <- read.table("old_id_map.tsv",sep="\t",header=T)
  new.data <- merge(data,old_ids,how="left",on="Protein")
  new.data <- new.data[-which(duplicated(new.data$Locus)),]
  fold.change <- new.data$dC
  names(fold.change) <- new.data$Locus
  fold.change <- sort(fold.change,decreasing = T)
  if (module)
  {
    gsea <- gseMKEGG(fold.change,organism="ctx")
  } else{
    gsea <- gseKEGG(fold.change,organism="ctx")
  }
  write.table(gsea@result,file.path(target.directory,ifelse(module,"kegg_gsea_module.tsv","kegg_gsea.tsv")),sep = "\t",col.names = T,row.names = F,quote=F)
}


diffExp.enrichment("DiffExp/Tables/hpt_vs_LL1210_early_log.tsv","DiffExp/KEGG/EL/",module = F)
diffExp.enrichment("DiffExp/Tables/hpt_vs_LL1210_early_log.tsv","DiffExp/KEGG/EL/",module = T)

diffExp.enrichment("DiffExp/Tables/hpt_vs_LL1210_mid_log.tsv","DiffExp/KEGG/ML/",module = F)
diffExp.enrichment("DiffExp/Tables/hpt_vs_LL1210_mid_log.tsv","DiffExp/KEGG/ML/",module = T)

diffExp.enrichment("DiffExp/Tables/hpt_vs_LL1210_late_log.tsv","DiffExp/KEGG/LL/",module = F)
diffExp.enrichment("DiffExp/Tables/hpt_vs_LL1210_late_log.tsv","DiffExp/KEGG/LL/",module = T)


diffExp.gsea("DiffExp/Tables/hpt_vs_LL1210_early_log.tsv","DiffExp/KEGG/EL/",module = F)
diffExp.gsea("DiffExp/Tables/hpt_vs_LL1210_early_log.tsv","DiffExp/KEGG/EL/",module = T)

diffExp.gsea("DiffExp/Tables/hpt_vs_LL1210_mid_log.tsv","DiffExp/KEGG/ML/",module = F)
diffExp.gsea("DiffExp/Tables/hpt_vs_LL1210_mid_log.tsv","DiffExp/KEGG/ML/",module = T)

diffExp.gsea("DiffExp/Tables/hpt_vs_LL1210_late_log.tsv","DiffExp/KEGG/LL/",module = F)
diffExp.gsea("DiffExp/Tables/hpt_vs_LL1210_late_log.tsv","DiffExp/KEGG/LL/",module = T)

diffCoExpress.enrichment("DiffCoExpress/differential_coexpression_strains.tsv",target.directory = "DiffCoExpress/KEGG/")
diffCoExpress.enrichment("DiffCoExpress/differential_coexpression_strains.tsv",target.directory = "DiffCoExpress/KEGG/",module=T)

diffCoExpress.gsea("DiffCoExpress/differential_coexpression_strains.tsv",target.directory = "DiffCoExpress/KEGG/",module = F)
diffCoExpress.gsea("DiffCoExpress/differential_coexpression_strains.tsv",target.directory = "DiffCoExpress/KEGG/",module = T)


coExpress.clust.enrichment("Clust/hpt/Clusters_Objects.tsv","Clust/hpt/KEGG/")
coExpress.clust.enrichment("Clust/hpt/Clusters_Objects.tsv","Clust/hpt/KEGG/",T)
coExpress.clust.enrichment("Clust/LL1210/Clusters_Objects.tsv","Clust/LL1210/KEGG/")
coExpress.clust.enrichment("Clust/LL1210/Clusters_Objects.tsv","Clust/LL1210/KEGG/",T)



coExpress.WGCNA.enrichment("WGCNA/Modules/hpt/hpt_modules.tsv","WGCNA/Modules/hpt/KEGG/")
coExpress.WGCNA.enrichment("WGCNA/Modules/hpt/hpt_modules.tsv","WGCNA/Modules/hpt/KEGG/",T)
coExpress.WGCNA.enrichment("WGCNA/Modules/LL1210/ll1210_modules.tsv","WGCNA/Modules/LL1210/KEGG/")
coExpress.WGCNA.enrichment("WGCNA/Modules/LL1210/ll1210_modules.tsv","WGCNA/Modules/LL1210/KEGG/",T)


