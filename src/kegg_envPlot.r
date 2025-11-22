#!/bin/env Rscript
# author: ph-u
# script: kegg_envPlot.r
# desc: explot coloring dataset for KEGG maps
# in: Rscript kegg_envPlot.r
# out: res/kegg_envPlot--*.csv
# arg: 0
# date: 20251013

source("p_src.r")
diffORF = read.csv(paste0(pT[2],"differentialORF--all.csv"), row.names = 1, header = T)

#d.diff = diffORF[which(diffORF[,2]>-1 & (diffORF[,1]!=diffORF[,2])),]
d.diff = diffORF[which(diffORF[,1]!=diffORF[,2]),]
r0 = data.frame(ORF=row.names(d.diff), color=ifelse(d.diff[,2]<0,cBp[4],ifelse(d.diff[,2]>0,cBp[14],cBp[16])))

write.table(r0, paste0(pT[1],"kegg_envPlot.csv"), sep = "\t", col.names = F, row.names = F, quote = F)
