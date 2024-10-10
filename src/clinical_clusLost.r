#!/bin/env Rscript
# author: ph-u
# script: clinical_clusLost.r
# desc: map Lost genes with STRING cluster
# in: Rscript clinical_clusLost.r
# out: NA
# arg: 0
# date: 20240823

source("p_src.r")
library(lattice)

lGene = read.csv(paste0(pT[1],"clinical_gMap.csv"), header = T, row.names = 1)
colnames(lGene) = read.table(text = gsub("_ASM","@",colnames(lGene)), sep = "@")[,1]
i.0 = mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==unique(mEta$sOurce)[grep("Cystic",unique(mEta$sOurce))])]
stringClus = read.table(paste0(pT[4],"string_kmeans_clusters--PA.tsv"), sep = "\t", header = T, comment.char = "", quote = "")

##### data trim: CF & within lost gene region #####
lGene = lGene[which(row.names(lGene) %in% stringClus$protein.name),which(colnames(lGene) %in% i.0)]
lG0 = lGene<1 # gene lost
print(stringClus[which(stringClus$protein.name %in% row.names(lG0)[which(rowSums(lG0)>(.7*ncol(lG0)))]),c(2,3,6)])

for(i in 1:nrow(lG0)){
    lG0[i,] = ifelse(lG0[i,],stringClus$cluster.number[which(stringClus$protein.name==row.names(lG0)[i])],0)
};rm(i)

pdf(paste0(pT[2],"clinical_clusLost--heat.pdf"), width = 10)
print(levelplot(t(lG0), col.regions = rep(c("#00000000",paste0(unique(stringClus$hex.color),"77", sep = "")), each = 4), xlab = "clinical isolates from individuals with Cystic Fibrosis", ylab = "PAO1 Gene Number (PA2125 -> PA2384)")) # cBp[c(3,5,6,7,4)]
invisible(dev.off())
