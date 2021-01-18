

phylo <- read.table("2020-07-22-phylo_prof_significant_pufs.tsv",sep="\t",header=T,stringsAsFactors = F)
phylo.sig <- phylo[which(phylo$Correlated.AICc - phylo$Independent.AICc <= -2),]

ortho <- read.table("Results_Jul22_3/Orthogroups/Orthogroups.tsv",sep="\t",header=T,stringsAsFactors = F)
ortho.ctherm <- ortho[,c("Orthogroup","cthermocellum_protein")]
phylo.sig.tmp <- merge(phylo.sig,ortho.ctherm,by.x="Protein",by.y="Orthogroup")
colnames(phylo.sig.tmp)[which(colnames(phylo.sig.tmp) == "cthermocellum_protein")] <- "Target_Protein_Ids" 
phylo.sig.tmp <- merge(phylo.sig.tmp,ortho.ctherm,by.x="PUF",by.y="Orthogroup")
colnames(phylo.sig.tmp)[which(colnames(phylo.sig.tmp) == "cthermocellum_protein")] <- "PUF_Protein_Ids" 

write.table(phylo.sig.tmp,"2020-07-22-phylo_prof_significant_pufs_with_prot_ids.tsv",sep="\t",col.names=T,row.names = F,quote=F)
