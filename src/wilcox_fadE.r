#!/bin/env Rscript
# author: ph-u
# script: wilcox_fadE.r
# desc: calculate Wilcox statistics for sequence windows, domain-specific
# in: Rscript wilcox_fadE.r
# out: res/wilcox_fadE.csv
# arg: 0
# date: 20240221

pT = paste0("../",c("data","res"),"/")
fadE1 = read.csv(paste0(pT[1],"01_PAO1_107_PA0506.csv"), header = T)
fadE2 = read.csv(paste0(pT[1],"01_PAO1_107_PA0508.csv"), header = T)
fadE1.0 = read.csv(paste0(pT[1],"stat_fadE1.csv"), header = T)
fadE2.0 = read.csv(paste0(pT[1],"stat_fadE2.csv"), header = T)

rEs0 = c("B","N","C1","C2")
rEs1 = c("Domain","fadE1.median","fadE2.median","Wilcox","p.val","p.adj")
rEs = as.data.frame(matrix(0,nr = length(rEs0), nc = length(rEs1)))
colnames(rEs) = rEs1; rEs$Domain = rEs0

for(i in 1:nrow(rEs)){
    e1 = fadE1$dNdS[which(fadE1.0[,which(colnames(fadE1.0)==rEs$Domain[i])]>0)]
    e2 = fadE2$dNdS[which(fadE2.0[,which(colnames(fadE2.0)==rEs$Domain[i])]>0)]
    e1 = e1[which(!is.na(e1))]
    e2 = e2[which(!is.na(e2))]
    w = wilcox.test(e1,e2)
    rEs[i,-1] = c(median(e1),median(e2),w$statistic,w$p.value,0)
};rm(i, w, e1, e2)

rEs$p.adj = p.adjust(rEs$p.val, "BH")

write.csv(rEs, paste0(pT[2],"wilcox_fadE.csv"), row.names = F, quote = F)
