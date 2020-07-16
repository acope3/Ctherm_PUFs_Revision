library(limma)
library(Biobase)
library(preprocessCore)
library(cowplot)
library(dplyr)
library(topGO)

createDFfromFit <- function(prots,fit)
{
  tmp <- data.frame(Protein=prots,Fold.Change=fit$coefficients[,1],P.value = fit$F.p.value,stringsAsFactors = F)
  return(tmp)
}

pufs <- read.table("../Data/all_pufs_dufs.txt",sep="",header=F,stringsAsFactors = F)
pufs <- pufs[,1]

prot.abund <- read.table("../Data/Expression/PUF_time_course_hpt_LLCtherm_min_2_prot_log2_Imputed_remove_crap.tsv",sep="\t",header=T,row.names=1,stringsAsFactors = F)
prot.abund.norm <- as.data.frame(normalize.quantiles(as.matrix(prot.abund)))
rownames(prot.abund.norm) <- rownames(prot.abund)
colnames(prot.abund.norm) <- colnames(prot.abund)
cols <- colnames(prot.abund)
sample.names <- unlist(strsplit(cols,split = "_rep[0-9]+"))
f <- factor(cols,levels=cols,labels=sample.names)
design <- model.matrix(~0+f)
colnames(design) <- unique(sample.names)
fit <- lmFit(prot.abund.norm,design)


pdf("Images/c_therm_limma_diff_expression_analysis.pdf")

cont.wt <- makeContrasts("LL1210_T1-hpt_T1",levels=design)
fit.el <- contrasts.fit(fit,cont.wt)
fit.el <- eBayes(fit.el)



tmp <- createDFfromFit(rownames(prot.abund),fit.el)
tmp["Color"] <- apply(tmp,MARGIN = 1,function(x){ifelse((x["Protein"] %in% pufs) && (as.numeric(x["P.value"]) < 0.05) && (abs(as.numeric(x["Fold.Change"])) > 2),"red","black")})

volcanoplot(fit.el,main=expression(Delta*"hpt vs. LL1210: Early-Log"),col=tmp$Color,pch=1,cex=1)
abline(v = -2,col="red",lty=2)
abline(v= 2,col="red",lty=2)
abline(h= -log10(0.05),col="red",lty=2)
legend(x="topleft",col = c("red"),pch = 1,legend=c("PUF"))

cont.wt <- makeContrasts("LL1210_T2-hpt_T2",levels=design)
fit.ml <- contrasts.fit(fit,cont.wt)
fit.ml <- eBayes(fit.ml)


tmp <- createDFfromFit(rownames(prot.abund),fit.ml)
tmp["Color"] <- apply(tmp,MARGIN = 1,function(x){ifelse((x["Protein"] %in% pufs) && (as.numeric(x["P.value"]) < 0.05) && (abs(as.numeric(x["Fold.Change"])) > 2),"red","black")})



volcanoplot(fit.ml,main=expression(Delta*"hpt vs. LL1210: Mid-Log"),col=tmp$Color,pch=1,cex=1)
abline(v = -2,col="red",lty=2)
abline(v= 2,col="red",lty=2)
abline(h= -log10(0.05),col="red",lty=2)
legend(x="topleft",col = c("red"),pch = 1,legend=c("PUF"))


cont.wt <- makeContrasts("LL1210_T4-hpt_T3",levels=design)
fit.ll <- contrasts.fit(fit,cont.wt)
fit.ll <- eBayes(fit.ll)


tmp <- createDFfromFit(rownames(prot.abund),fit.ll)
tmp["Color"] <- apply(tmp,MARGIN = 1,function(x){ifelse((x["Protein"] %in% pufs) && (as.numeric(x["P.value"]) < 0.05) && (abs(as.numeric(x["Fold.Change"])) > 2),"red","black")})


volcanoplot(fit.ll,main=expression(Delta*"hpt vs. LL1210: Late-Log"),col=tmp$Color,pch=1,cex=1)
abline(v = -2,col="red",lty=2)
abline(v= 2,col="red",lty=2)
abline(h= -log10(0.05),col="red",lty=2)
legend(x="topleft",col = c("red"),pch = 1,legend=c("PUF"))


dev.off()

fit.el.table <- topTable(fit.el,number = nrow(prot.abund))
fit.ml.table <- topTable(fit.ml,number = nrow(prot.abund))
fit.ll.table <- topTable(fit.ll,number = nrow(prot.abund))


fit.el.table["PUF"] <- ifelse(rownames(fit.el.table) %in% pufs,"Yes","No")
fit.ml.table["PUF"] <- ifelse(rownames(fit.ml.table) %in% pufs,"Yes","No")
fit.ll.table["PUF"] <- ifelse(rownames(fit.ll.table) %in% pufs,"Yes","No")

write.table(fit.el.table,"Tables/hpt_vs_LL1210_early_log.tsv",row.names = T,col.names = T,quote=F,sep="\t")
write.table(fit.ml.table,"Tables/hpt_vs_LL1210_mid_log.tsv",row.names = T,col.names = T,quote=F,sep="\t")
write.table(fit.ll.table,"Tables/hpt_vs_LL1210_late_log.tsv",row.names = T,col.names = T,quote=F,sep="\t")



#t.prot.abund <- as.data.frame(t(prot.abund))

# cols <- rownames(t.prot.abund)
# sample.names <- unlist(strsplit(cols,split = "_rep[0-9]+"))
# t.prot.abund["Sample.names"] <- sample.names
# 
# x <- t.prot.abund %>% group_by(Sample.names) %>% summarise_all(mean)
# x <- as.data.frame(t(x))
# colnames(x) <- as.character(unlist(unname(x[1,])))
# x <- x[!(rownames(x) %in% c("Sample.names")),]
# x[] <- lapply(x, function(i) as.numeric(as.character(i)))
# 
# colramp = colorRampPalette(c(3,"white",2))(ncol(x))
# plot(density(x[,1]),col=colramp[1],lwd=3,main="Density plot for Mean Protein Abundance\nBy Sample Type")
# for(i in 2:ncol(x))
# {
#   lines(density(x[,i]),lwd=3,col=colramp[i])
# }
# legend(x = "topright",legend = colnames(x),col = colramp,lty=1)
# 
# norm.x <- normalize.quantiles(as.matrix(x))
# plot(density(norm.x[,1]),col=colramp[1],lwd=3,main="Density plot for Normalized Mean Protein Abundance\nBy Sample Type")
# for(i in 2:ncol(norm.x))
# {
#   lines(density(norm.x[,i]),lwd=3,col=colramp[i])
# }
# legend(x = "topright",legend = colnames(x),col = colramp,lty=1)




