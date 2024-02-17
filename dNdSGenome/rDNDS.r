#!/bin/env Rscript
# author: ph-u
# script: rDNDS.r
# desc: calculate & plot residue-wise dN/dS contribution of the given gene
# in: Rscript rDNDS.r [srcFile] [prefix] [gene] [locus_Tag]
# out: res/02_[prefix]_[gene]_[locus_Tag].{csv,pdf}
# arg: 4
# date: 20240119

argv=(commandArgs(T))
# argv = c("../data/01_eColiK12_b0462.csv", "eColiK12", "b0462", "b0462")
# argv=c("../data/01_PAO1_107_PGD134012.csv","PAO1_107","PGD134012","PA0001")
a = read.csv(argv[1], header = T)
a = a[which(!is.na(a$dNdS)),]

# ranged from 0 (all dS) to inf (all dN)
aMino = data.frame(ntPos=seq(1, max(c(a$start,a$end)), 3), mEdian=NA, nRef=NA, p.val=NA)

##### max number of available dNdS values #####
sqMax = 1:nrow(aMino)
sqMax[(ceiling(nrow(aMino)/2)+!(nrow(aMino)%%2)):nrow(aMino)] = ceiling(nrow(aMino)/2):1
sqMax[sqMax>(a$end[1]-a$start[1])] = (a$end[1]-a$start[1])
sqMax = sqMax*length(list.files(".", ".ndb"))

##### Protein residue essentiality potential #####
for(i in 1:nrow(aMino)){
    i0 = which(a$start <= aMino$ntPos[i] & a$end >= aMino$ntPos[i])
    aMino$mEdian[i] = median(a$dNdS[i0])
    aMino$nRef[i] = length(i0)
    aMino$p.val = wilcox.test(a$dNdS[i0], a$dNdS)$p.value
};rm(i, i0)
aMino$p.adj = p.adjust(aMino$p.val, method = "BH")
aMino$mEdian[which(aMino$p.adj > .05)] = median(a$dNdS)
aMino$cOl = gray(1-(.1+.9*aMino$nRef/sqMax))

pdf(paste0("../res/02_",argv[2],"_",argv[3],"_",argv[4],".pdf"), width = 21)
plot(1:nrow(aMino), aMino$mEdian, col=aMino$cOl, type="p", pch=19, xlab = "Residue position", ylab = "Reconstructed Median dN/dS", main = paste(argv[3],"/",argv[4]), ylim = c(0, max(1,aMino$mEdian)))
lines(x=1:nrow(aMino), y=aMino$mEdian, col="#000000ff", lty=1)
abline(h=1, col="#00ccff77", lty=2, lwd=3)
invisible(dev.off())

write.csv(aMino, paste0("../res/02_",argv[2],"_",argv[3],"_",argv[4],".csv"), row.names = F, quote = F)
