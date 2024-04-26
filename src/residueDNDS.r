#!/bin/env Rscript
# author: ph-u
# script: residueDNDS.r
# desc: Plot "--reCon.csv" mean dN/dS residue-wise resolution
# in: Rscript residueDNDS.r [gene]
# out: res/[gene]_reCon.pdf
# arg: 1
# date: 20240406, 20240424

argv = (commandArgs(T))
source("p_src.r")
pRe = c("Cystic fibrosis", "Environmental", "Other infections", "Skin", "Unknown")
stP = 14

for(i in 1:length(pRe)){
    a = read.csv(paste0(pT[1], "PAO1_107_",argv[1],"_", gsub(" ","-",pRe[i]), "--reCon.csv"), header = T)
    if(i==1){p0 = a[,c("codon", "aaRes", "dNdS.mean")]; colnames(p0)[3] = gsub(" ","-",pRe[i])}else{p0[,gsub(" ","-",pRe[i])] = a$dNdS.mean}
};rm(i, a)
#pdf(paste0(pT[2], argv[1],"_reCon.pdf"), width = 50, height = 7)
pdf(paste0(pT[2], argv[1],"_reCon.pdf"), width = 7, height = 9)
par(mar = c(4,4,1,1)+.1, mfrow = c(2,1), cex.axis = 1.2)
matplot(x = row.names(p0), y = p0[,-(1:2)], type = "b", pch = stP:(stP-1+length(pRe)), col = cBp[1:length(pRe)], xlab = paste0(argv[1]," Residue position"), ylab = "Single codon resolution mean dN/dS", cex = .7)
abline(h=1, col = "#000000aa", lty = 3, lwd = 2)
#axis(1, at = 1:nrow(p0), labels = paste0(p0$codon, "\n[", p0$aaRes, "]"))
#axis(1, at = 1:nrow(p0), labels = p0$aaRes)
plot.new()
legend("top", legend = pRe, col = cBp[1:length(pRe)], pch = stP:(stP-1+length(pRe)), cex = 1.5, bg = "#ffffffff")
invisible(dev.off())
