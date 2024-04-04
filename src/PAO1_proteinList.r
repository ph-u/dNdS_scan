#!/bin/env Rscript
# author: ph-u
# script: PAO1_proteinList.r
# desc: plot cumulative protein list unique entries
# in: Rscript PAO1_proteinList.r
# out: res/PAO1_proteinList.pdf
# arg: 0
# date: 20240317

pT = c(paste0("../",c("data","res"),"/"), "/media/pokman/Transcend/PAO1_dNdS/data/")
f = list.files(pT[3],"^01_")
#a = read.csv(paste0(pT[1],"PAO1_proteins.txt"), header = F)
a0 = read.csv("../p_pdf/PAO1_proteomics/PAO1_proteins.csv", header = F)

rEs = data.frame(no=0:ncol(a0), prop=0)
for(i in 2:nrow(rEs)){
    r0 = unique(unname(unlist(strsplit(unname(unlist(a0[,1:rEs$no[i]])),";"))))
    rEs$prop[i] = round(sum(substr(r0,1,2)=="PA")/length(f)*100)
};rm(i, r0)

pdf(paste0(pT[2],"PAO1_proteinList.pdf"))
par(mar = c(4,4,0,0)+.1)
plot(rEs[,1], rEs[,2], "b", xlab = "Number of selected publications", ylab = "Cumulative Percentage of PAO1 genes recorded in proteomics (%)", ylim = c(0,100), pch = 3)
invisible(dev.off())
