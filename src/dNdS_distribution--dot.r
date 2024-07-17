#!/bin/env Rscript
# author: ph-u
# script: dNdS_distribution--dot.r
# desc: plot mean dN/dS distribution, CF vs environmental
# in: Rscript dNdS_distribution--dot.r
# out: res/dNdS_distribution--dot*
# arg: 0
# date: 20240519

source("p_src.r");source("p_metabolism_PAO1.r")
srCount = table(mEta$sOurce)
a = read.csv(paste0(pT[1], "dNdS_distribution.csv"), header = T)
a$mean = ifelse(is.infinite(a$mean), ifelse(is.infinite(a$Q2),ifelse(is.infinite(a$Q3),ifelse(is.infinite(a$Q1),ifelse(is.infinite(a$max),max(a$mean[is.finite(a$mean)], na.rm = T)+2,a$max),a$Q1),a$Q3),a$Q2), a$mean) # median replace inf mean
a0 = data.frame(gene = unique(a$gene), CF = a$mean[which(a$sOurce=="Cystic fibrosis")], env = a$mean[which(a$sOurce=="Environmental")], CF.NA = round(a$NAs[which(a$sOurce=="Cystic fibrosis")]/srCount[which(names(srCount)=="Cystic fibrosis")], 2), env.NA = round(a$NAs[which(a$sOurce=="Environmental")]/srCount[which(names(srCount)=="Environmental")], 2))

a1 = a0[which(a0$gene %in% gBioc$protList),]
a1$cOl = gBioc$cOl[match(a1$gene, gBioc$protList)]
pdf(paste0(pT[2], "dNdS_distribution--dot.pdf"), width = 8, height = 10)
par(mar = c(4,4,1,0)+.1, mfrow = c(3,2))
for(i in 1:2){
    plot(x = a1$CF, y = a1$env, col = a1$cOl, pch = 19, xlab = "Mean dN/dS in CF", ylab = "Mean dN/dS in environment", cex = .2)
    if(i==2){text(x = a0$CF, y = a0$env, labels = a0$gene, cex = .5, col = cBp[length(cBp)])}
    abline(a=0, b=1, h=1, v=1, lty = 2, lwd = 1)
}
plot(x = a1$CF, y = a1$env, col = a1$cOl, pch = 19, xlab = "Mean dN/dS in CF", ylab = "Mean dN/dS in environment", xlim = c(0,1), ylim = c(0,1), cex = .2)
#text(x = a0$CF, y = a0$env, labels = a0$gene, cex = .1, col = cBp[length(cBp)])
abline(a=0, b=1, h=1, v=1, lty = 2, lwd = 1)
plot.new()
mB = unique(gBioc[,-1])
legend("center", pch = 19, legend = mB[,1], col = mB[,2], ncol = 2, cex = .5)

for(i in 1:2){
    plot(x = log10(a1$CF+.00001), y = log10(a1$env+.00001), col = a1$cOl, pch = 19, xlab = "log10 Mean dN/dS in CF", ylab = "log10 Mean dN/dS in environment", cex = .2)
    if(i==2){text(x = log10(a0$CF+.00001), y = log10(a0$env+.00001), labels = a0$gene, cex = .5, col = cBp[length(cBp)])}
    abline(a=0, b=1, h=0, v=0, lty = 2, lwd = 1)
};rm(i)
invisible(dev.off())

a0$wtinRange = apply(a0[,2:3], 1, diff)<.5
fEat = read.table("../raw/features.txt", sep = "\t", header = T, quote = "")
fEat.0 = fEat[match(a0$gene, fEat$Locus.Tag),c("Gene.Name", "Product.Description")]
write.table(cbind(a0, fEat.0), paste0(pT[1],"dNdS_distribution--dotFull.csv"), sep = "\t", row.names = F, col.names = T, quote = F)

write.table(a0[which(a0$wtinRange==F),], paste0(pT[2],"dNdS_distribution--dot.csv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(a0[which(a0$CF>2 | a0$env>2),], paste0(pT[2],"dNdS_distribution--dot-highPos.csv"), sep = "\t", row.names = F, col.names = T, quote = F)


