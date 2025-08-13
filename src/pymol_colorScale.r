#!/bin/env Rscript
# author: ph-u
# script: pymol_colorScale.r
# desc: rescale dN/dS values (PA1430, PA1373)
# in: Rscript pymol_colorScale.r
# out: data/[PAnum]_{Cystic-fibrosis,Environmental}.csv
# arg: 0
# date: 20250717

PAnum = paste0("PA",c(5086:5088)) # c(1430,2373)
sRc = c("Cystic-fibrosis","Environmental")

for(i in 1:length(PAnum)){
  d.c = read.csv(paste0("../data/PAO1_107_",PAnum[i],"_",sRc[1],"--reCon.csv"), header = T)
  d.e = read.csv(paste0("../data/PAO1_107_",PAnum[i],"_",sRc[2],"--reCon.csv"), header = T)

  d.c$dNdS.mean[nrow(d.c)-1] = d.e$dNdS.mean[nrow(d.e)-1] = max(c(d.c$dNdS.mean,d.e$dNdS.mean))

  write.csv(d.c, paste0("../data/",PAnum[i],"_",sRc[1],".csv"), row.names = F, quote = F)
  write.csv(d.e, paste0("../data/",PAnum[i],"_",sRc[2],".csv"), row.names = F, quote = F)
}

