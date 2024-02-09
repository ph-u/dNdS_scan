#!/bin/env Rscript
# author: ph-u
# script: aaAlign.r
# desc: align amino acid sequences
# in: Rscript aaAlign.r
# out: data/*.aln
# arg: 0
# date: 20240204

library(ape)
library(msa)
pT = c("../data/")

f = list.files(pT[1],"_dbAA.fa")
for(i in 1:length(f)){
    a = msa(readAAStringSet(paste0(pT[1],f[i])), "ClustalOmega")
#    print(a, show = "complete")
    writeXStringSet(unmasked(a), file=paste0(pT[1],gsub("_dbAA", "-aaAlign", f[i])))
};rm(i)
