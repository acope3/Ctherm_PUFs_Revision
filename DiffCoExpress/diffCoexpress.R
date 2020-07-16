library(DCGL)


hpt <- read.table("Cleaned_data/hpt_cleaned_imputed_data_transposed.tsv",sep="\t",header=T,stringsAsFactors = F)
ll <- read.table("Cleaned_data/ll1210_cleaned_imputed_data_transposed.tsv",sep="\t",header = T,stringsAsFactors = F)

expr <- merge(hpt,ll,on="Protein",how="inner")

expr.A <- expr[,2:13]
rownames(expr.A) <- expr$Protein

expr.B <- expr[,14:29]
rownames(expr.B) <- expr$Protein

DCp.res <- DCp(expr.A,expr.B,cutoff = 0.25,N=100)


diff.coex <- WGCNA(expr.A,expr.B,variant = "DCp")
