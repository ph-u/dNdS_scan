#!/bin/env Rscript
# author: ph-u
# script: dNdS_distribution.r
# desc: get dN/dS distribution of every gene in PAO1 relative to clinical isolates, source dependent
# in: Rscript dNdS_distribution.r
# out: data/dNdS_distribution.csv
# arg: 0
# date: 20240517

source("p_src.r")
f.dbSum = list.files(pT[3], "--dbSum")
f.dbSum = data.frame(fNam = f.dbSum, gEne = read.table(text=read.table(text=f.dbSum, sep = "_")[,3], sep = "-")[,1])

cat("Summarizing gene dN/dS distributions by source:", date(), "\n")
for(i in 1:nrow(f.dbSum)){
    cat("Summarizing (",i,") gene", f.dbSum$gEne[i], ":", date(), "     \r")
    i0 = read.csv(paste0(pT[3], f.dbSum[i,1]), header = T)
    i0$clinical = read.table(text=gsub("_ASM", "@", i0$clinical), sep = "@")[,1]
    i0$sRc = mEta$sOurce[match(i0$clinical, mEta$assemblyInfo.genbankAssmAccession)]
    i1 = tapply(i0$dNdS, i0$sRc, summary)
    r0.c = c("gene", "sOurce", "min", "Q1", "Q2", "mean", "Q3", "max", "NAs")
    r0 = as.data.frame(matrix(nr = length(i1), nc = length(r0.c)))
    colnames(r0) = r0.c
    r0[,1] = f.dbSum$gEne[i]
    r0[,2] = names(i1)
    for(i2 in 1:length(i1)){ if(length(i1[[i2]])==2){ r0[i2, -c(1:2)] = c(rep(NA, ncol(r0)-3), as.numeric(i1[[i2]][2])) }else if(length(i1[[i2]])<7){ r0[i2, -c(1:2)] = c(i1[[i2]],0) }else{r0[i2, -c(1:2)] = i1[[i2]]} }
    if(i>1){
        write.csv(rbind(read.csv(paste0(pT[1],"dNdS_distribution.csv"), header = T), r0), paste0(pT[1],"dNdS_distribution.csv"), row.names = F, quote = F)
    }else{
        write.csv(r0, paste0(pT[1],"dNdS_distribution.csv"), row.names = F, quote = F)
    }
};rm(i, i0, i1, i2, r0);cat("\nDone summary by sampling source: ", date(), "\n")

##### export CF & env, reorder dN/dS in decending order #####
#r0 = read.csv("../data/dNdS_distribution.csv", header = T)
#for(i in c("Cystic fibrosis", "Environmental")){
#    r1 = r0[which(r0$sOurce==i),]
#    write.csv(r1[rev(order(r1$mean)),], paste0(pT[1], "dNdS_distribution_",substr(i,1,3),".csv"), row.names = F, quote = F)
#};rm(i)
