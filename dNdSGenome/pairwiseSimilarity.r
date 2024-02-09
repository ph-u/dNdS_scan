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
#library(protr) # https://cran.r-project.org/web/packages/protr/vignettes/protr.html

argv = list.files("../data","_db.fa")[as.numeric(argv)]
library(stringdist)
library(ape)
source("src_dNdS.r")

a = as.character(read.FASTA(paste0("../data/",argv[1])))
aSpace = numeric(length(a))
for(i in 1:length(a)){
    aSp0 = which(a[[i]]=="-")
    aSpace[i] = length(aSp0)
    a[[i]][aSp0] = ""
    a[[i]] = strsplit(paste0(a[[i]], collapse = ""), "")[[1]]
};rm(i, aSp0)

p = read.csv(paste0("../data/",gsub("_db.fa",".csv",argv[1])), header = T)
a0 = vector(mode="list", length = length(a))
p0 = as.data.frame(matrix(0, nr = length(a), nc = 2))
names(a0) = names(a)
for(i in 1:length(a)){
    i0 = as.numeric(p[grep(strsplit(names(a0)[i], ";")[[1]][1], p$clinical)[1],c("sample_start", "len")]) - c(0,aSpace[i])
    # trim nucleotide overhang to match ORF
    i0 = c(i0[1], i0[1] + (4-ifelse(i0[1]%%3==0,3,i0[1]%%3))%%3, i0[2] - (i0[1]+i0[2]-1)%%3)
    a0[[i]] = nt2prot(paste0(a[[i]][(i0[2]-i0[1]+1):i0[3]], collapse = ""))
    p0[i,] = c(ceiling(i0[2]/3),nchar(a0[[i]]))
};rm(i)
write.FASTA(as.AAbin(a0), paste0("../data/",gsub("_db","_dbAA",argv[1])))
#parSeqSim(a0)

x = c("seq1","seq2","simRatio")
a1 = as.data.frame(matrix(nr = choose(length(a),2), nc = length(x)))
colnames(a1) = x; rm(x)
date(); i = 1; for(i0 in 1:(length(a)-1)){ for(i1 in (i0+1):length(a)){
    cat("Pairwise compare: ", i, " / ",nrow(a1)," - ",date(),"\n") #"\r")
    if(all(p0[i0,] == p0[i1,])){
        s1 = a0[[i0]]
	s2 = a0[[i1]]
    }else{
        s0 = c(max(p0[i0,1],p0[i1,1]), min(p0[i0,2],p0[i1,2]))
        s1 = s0[1]-p0[i0,1]+1
        s2 = s0[1]-p0[i1,1]+1
        s1 = substr(a0[[i0]],s1,s1+s0[2])
	s2 = substr(a0[[i1]],s2,s2+s0[2])
	rm(s0)
    }
    a1[i,] = c(i0,i1,stringsim(s1,s2))
    i = i+1
}};rm(i, i0, i1, s1, s2)
a1$rAnge = round(a1$simRatio,2)
write.csv(a1,paste0("../data/",gsub(".fa$",".csv",argv[1])), quote = F, row.names = F)
