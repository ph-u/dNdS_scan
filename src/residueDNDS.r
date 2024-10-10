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
pRe = c("Cystic fibrosis", "Environmental", "Other infections", "Urban", "Unknown")
#pRe = c("Cystic fibrosis", "Environmental")
stP = 14

for(i in 1:length(pRe)){
    a = read.csv(paste0(pT[1], "PAO1_107_",argv[1],"_", gsub(" ","-",pRe[i]), "--reCon.csv"), header = T)
    if(i==1){
        p0 = a[,c("codon", "aaRes", "dNdS.1Q")]
        p1 = a[,c("codon", "aaRes", "dNdS.mean")]
        p2 = a[,c("codon", "aaRes", "dNdS.2Q")]
        p3 = a[,c("codon", "aaRes", "dNdS.3Q")]
        colnames(p0)[3] = colnames(p1)[3] = colnames(p2)[3] = colnames(p3)[3] = gsub(" ","-",pRe[i])
    }else{
        p0[,gsub(" ","-",pRe[i])] = a$dNdS.1Q
        p1[,gsub(" ","-",pRe[i])] = a$dNdS.mean
        p2[,gsub(" ","-",pRe[i])] = a$dNdS.2Q
        p3[,gsub(" ","-",pRe[i])] = a$dNdS.3Q
    }
};rm(i, a)
#pdf(paste0(pT[2], argv[1],"_reCon.pdf"), width = 50, height = 7)
pdf(paste0(pT[2], argv[1],"_reCon.pdf"), width = 7, height = 14)
par(mar = c(6,6,1,1)+.1, mfrow = c(5,1), cex.axis = 3, cex.lab = 3)

matplot(x = row.names(p0), y = p0[,-(1:2)], type = "b", pch = stP:(stP-1+length(pRe)), col = cBp[1:length(pRe)], xlab = paste0(argv[1]," Residue position"), ylab = "dN/dS, 1Q", cex = .7)
abline(h=1, col = "#000000aa", lty = 3, lwd = 2)

matplot(x = row.names(p1), y = p1[,-(1:2)], type = "b", pch = stP:(stP-1+length(pRe)), col = cBp[1:length(pRe)], xlab = paste0(argv[1]," Residue position"), ylab = "dN/dS, mean", cex = .7)
abline(h=1, col = "#000000aa", lty = 3, lwd = 2)

matplot(x = row.names(p2), y = p2[,-(1:2)], type = "b", pch = stP:(stP-1+length(pRe)), col = cBp[1:length(pRe)], xlab = paste0(argv[1]," Residue position"), ylab = "dN/dS, 2Q", cex = .7)
abline(h=1, col = "#000000aa", lty = 3, lwd = 2)

matplot(x = row.names(p3), y = p3[,-(1:2)], type = "b", pch = stP:(stP-1+length(pRe)), col = cBp[1:length(pRe)], xlab = paste0(argv[1]," Residue position"), ylab = "dN/dS, 3Q", cex = .7)
abline(h=1, col = "#000000aa", lty = 3, lwd = 2)

#axis(1, at = 1:nrow(p0), labels = paste0(p0$codon, "\n[", p0$aaRes, "]"))
#axis(1, at = 1:nrow(p0), labels = p0$aaRes)

plot.new()
legend("top", legend = pRe, col = cBp[1:length(pRe)], pch = stP:(stP-1+length(pRe)), cex = 1.5, bg = "#ffffffff")
invisible(dev.off())
