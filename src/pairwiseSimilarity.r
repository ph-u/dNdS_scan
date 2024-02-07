#!/bin/env Rscript
# author: ph-u
# script: pairwiseSimilarity.r
# desc: pairwise similarity score calculation for FASTA
# in: Rscript pairwiseSimilarity.r
# out: NA
# arg: 0
# date: 20240122

argv = (commandArgs(T))
# argv = paste0("../data/",c("00_eColiK12_b0462_db.fa","00_eColiK12_b0463_db.fa","00_eColiK12_b0762_db.fa","00_eColiK12_b3035_db.fa"))
library(stringdist)
#library(protr) # https://cran.r-project.org/web/packages/protr/vignettes/protr.html
library(ape)
p=getwd()
setwd("../../1_04_wholeGenomeDNDS/dNdSGenome")
source("src_dNdS.r")
setwd(p)

a = as.character(read.FASTA(paste0("../data/",argv[1],".fa")))
a0 = vector(mode="list", length = length(a))
names(a0) = names(a)
for(i in 1:length(a)){a0[[i]] = nt2prot(paste0(a[[i]], collapse = ""))};rm(i)
#parSeqSim(a0)

x = c("seq1","seq2","simRatio")
a1 = as.data.frame(matrix(nr = choose(length(a),2), nc = length(x)))
colnames(a1) = x; rm(x)
date(); i = 1; for(i0 in 1:(length(a)-1)){ for(i1 in (i0+1):length(a)){
    cat(i, " / ",nrow(a1)," - ",date(),"       \r")
    a1[i,] = c(i0,i1,stringsim(a0[[i0]],a0[[i1]]))
    i = i+1
}};rm(i0,i1)
a1$rAnge = round(a1$simRatio,2)
write.csv(a1,paste0("../data/",argv[1],".csv"), quote = F, row.names = F)
