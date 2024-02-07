#!/bin/env Rscript
# author: ph-u
# script: nFlipSeq.r
# desc: convert minus sequence into comparable complementary sequence
# in: Rscript nFlipSeq.r [atcg]
# out: stdout message
# arg: 1
# date: 20240105

argv=(commandArgs(T))
source("src_dNdS.r")

a = strsplit(argv[1], "")[[1]]
for(i in 1:length(a)){a[i] = nFlip(a[i],"-")}
cat(toupper(paste0(a, collapse = "")))
