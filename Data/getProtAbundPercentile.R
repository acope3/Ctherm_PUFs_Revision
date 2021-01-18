library(dplyr)

calculateMeanOverRep <- function(prot.abund)
{
  t.prot.abund <- as.data.frame(t(prot.abund))
  cols <- rownames(t.prot.abund)
  sample.names <- unlist(strsplit(cols,split = "_rep[0-9]+"))
  t.prot.abund["Sample.names"] <- sample.names
  
  x <- t.prot.abund %>% group_by(Sample.names) %>% summarise_all(mean)
  x <- as.data.frame(t(x))
  colnames(x) <- as.character(unlist(unname(x[1,])))
  x <- x[!(rownames(x) %in% c("Sample.names")),]
  x[] <- lapply(x, function(i) as.numeric(as.character(i)))
  return(x)
}

getMostAbundant <- function(df,row.index,percent=0.9)
{
  cutoff <- quantile(df[,row.index],probs = percent)
  df.abund <- df[which(df[,row.index] > cutoff),]
  return(rownames(df.abund))
}


hpt <- read.table("Expression/PUF_time_course_hpt_min_2_prot_log2_Imputed_remove_crap.tsv",sep="\t",header=T,stringsAsFactors = F,row.names=1)
ll <- read.table("Expression/PUF_time_course_LL1210_min_2_prot_log2_Imputed_remove_crap.tsv",sep="\t",header=T,stringsAsFactors = F,row.names = 1)

hpt.mean <- calculateMeanOverRep(hpt)
ll.mean <- calculateMeanOverRep(ll)

pufs <- read.table("measured_pufs.txt",stringsAsFactors = F)
pufs <- pufs[,1]

abundant.pufs <- vector(mode = "list",length=7)

for (i in 1:length(colnames(hpt.mean)))
{
  most.abund <- getMostAbundant(hpt.mean,i)
  abund.pufs <- pufs[which(pufs %in% most.abund)]
  abundant.pufs[[i]] <- abund.pufs
}

hpt.pufs <- Reduce(intersect,abundant.pufs[1:3])

for (i in 1:length(colnames(ll.mean)))
{
  most.abund <- getMostAbundant(ll.mean,i)
  abund.pufs <- pufs[which(pufs %in% most.abund)]
  abundant.pufs[[i+3]] <- abund.pufs
}

ll.pufs <- Reduce(intersect,abundant.pufs[4:7])

all.pufs <- Reduce(intersect,abundant.pufs)