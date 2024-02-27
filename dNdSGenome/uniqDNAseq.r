#!/bin/env Rscript
# author: ph-u
# script: uniqDNAseq.r
# desc: extract unique DNA sequences from database FASTA file
# in: Rscript uniqDNAseq.r [strain] [gene]
# out: data/00_[strain]_[gene]_db.fa
# arg: 2
# date: 20240227

argv = (commandArgs(T))
library(ape)
pT = "../data/"
nAm = paste0(pT,"00_",argv[1],"_",argv[2],"_db.fa")
a = as.character(read.FASTA(nAm))

aNam = unique(names(a))
a0 = vector(mode="list", length = length(aNam))
names(a0) = aNam
for(i in 1:length(aNam)){
    aTmp = which(names(a)==names(a0)[i])
    a0[[i]] = a[[aTmp[length(aTmp)]]] # extract latest extracted sequence
};rm(i)

write.FASTA(as.DNAbin(a0), nAm)
