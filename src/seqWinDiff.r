#!/bin/env Rscript
# author: ph-u
# script: seqWinDiff.r
# desc: Sequence windows variation by sources
# in: Rscript seqWinDiff.r
# out: NA
# arg: 1
# date: 20240425

argv = (commandArgs(T))
#argv = "PA2185"

source("p_src.r")
f = list.files(pT[3],argv[1])
rDNDS = read.csv(paste0(pT[3],f[grep("rDNDS",f)]), header = T)
rDNDS$sRc = mEta$sOurce[match(read.table(text = gsub("_ASM","@",rDNDS$clinical),sep = "@")[,1], mEta$assemblyInfo.genbankAssmAccession)]
print(pairwise.wilcox.test(rDNDS$dNdS, rDNDS$sRc, p.adjust = "bonf"))
