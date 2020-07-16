library(WGCNA)
library(Biobase)
library(preprocessCore)

## Based on WGCNA tutorials 

reformatMatrix <- function(df)
{
  datExpr0 <- as.data.frame(t(df))
  colnames(datExpr0) <- rownames(df)
  rownames(datExpr0) <- colnames(df)
  return(datExpr0)
}

## Quantile normalization
normalization <- function(df)
{
  df.norm <- as.data.frame(normalize.quantiles(as.matrix(df)))
  rownames(df.norm) <- rownames(df)
  colnames(df.norm) <- colnames(df)
  return(df.norm)
}


## change this file to switch between hpt and LL1210
hpt <- read.table("../Data/Expression/PUF_time_course_LL1210_min_2_prot_log2_remove_crap.tsv",header=T,row.names=1,sep="\t",stringsAsFactors = F)
hpt.imp <- read.table("../Data/Expression/PUF_time_course_LL1210_min_2_prot_log2_Imputed_remove_crap.tsv",header=T,row.names=1,sep="\t",stringsAsFactors = F)
hpt[hpt == 0] <- NA

hpt <- normalization(hpt)
hpt.imp <- normalization(hpt.imp)


hpt <- reformatMatrix(hpt)
hpt.imp <- reformatMatrix(hpt.imp)


gsg <- goodSamplesGenes(hpt, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(hpt)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(hpts)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data
  hpt <- hpt[gsg$goodSamples,gsg$goodGenes]
  hpt.imp <- hpt.imp[gsg$goodSamples,gsg$goodGenes]
}

sampleTree <- hclust(dist(hpt), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "Plots/sampleClustering_LL1210_not_imputed.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers in LL1210 non-imputed data", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()

sampleTree <- hclust(dist(hpt.imp), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "Plots/sampleClustering_LL1210_imputed.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers in LL1210 imputed data", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


write.table(hpt,"Cleaned_data/ll1210_cleaned_data.tsv",sep="\t",col.names = T,row.names = T,quote=F)
write.table(hpt.imp,"Cleaned_data/ll1210_cleaned_imputed_data.tsv",sep="\t",col.names = T,row.names = T,quote=F)


write.table(as.data.frame(t(hpt)),"Cleaned_data/ll11210_cleaned_data_transposed.tsv",sep="\t",col.names = T,row.names = T,quote=F)
write.table(as.data.frame(t(hpt.imp)),"Cleaned_data/ll1210_cleaned_imputed_data_transposed.tsv",sep="\t",col.names = T,row.names = T,quote=F)
