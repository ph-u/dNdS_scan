#!/bin/env Rscript
# author: ph-u
# script: density_dbSum.r
# desc: density plot on the whole-gene dN/dS ratio
# in: Rscript density_dbSum.r
# out: res/density_dbSum.jpeg
# arg: 0
# date: 20251012

source("p_src.r")
diffORF = read.csv(paste0(pT[2],"differentialORF--all.csv"), row.names = 1, header = T)

##### Collect dN/dS summary ##### (2 mins)
#r0.c = paste0("Q",1:3,".",rep(c("CF","environmental"), each = 3))
#r0 = as.data.frame(matrix(0, nrow = nrow(diffORF), ncol = length(r0.c)))
#colnames(r0) = r0.c
i1 = c("Cystic", "Envir")
#sEq = seq(1,ncol(r0),3)
cat(date(),": Calculating statistical summaries of ORFs\n")
d.cf = d.en = c()
for(i in 1:nrow(diffORF)){
  cat(date(),":",i,"/",nrow(diffORF),"(",round(i/nrow(diffORF)*100,2),"% )     \r")
  i0 = read.csv(paste0(pT[3],sub("_db.fa","--dbSum.csv",sub("00_","",f$fNam[which(f$gEne==row.names(diffORF)[i])]))), header = T)
  i0$clinical = read.table(text = sub("_ASM","@",i0$clinical), sep = "@")[,1]
  i0$dNdS = abs(i0$dNdS)
  d.cf = c(d.cf, i0$dNdS[which(i0$clinical %in% mEta$assemblyInfo.genbankAssmAccession[grep(i1[1], mEta$sOurce)])])
  d.en = c(d.en, i0$dNdS[which(i0$clinical %in% mEta$assemblyInfo.genbankAssmAccession[grep(i1[2], mEta$sOurce)])])
#  for(i2 in 1:length(i1)){
#    r0[i,sEq[i2]:(sEq[i2]+2)] = summary(i0$dNdS[which(i0$clinical %in% mEta$assemblyInfo.genbankAssmAccession[grep(i1[i2], mEta$sOurce)])], na.rm = T)[c("1st Qu.","Median","3rd Qu.")]
#  };rm(i2)
};rm(i,i0)
cat(date(),": Done statistical summaries of ORFs\n")

##### Plot IQR plot #####
#lIm = 3
#r1 = log10(r0+1/(10^lIm))
#i0 = unique(unlist(r0))[rev(order(unique(unlist(r0))))]
#i0 = log10(i0[!is.na(i0) & is.finite(i0)][1])+(lIm-1)
#for(i in 1:ncol(r1)){if(any(is.infinite(r1[,i]))){r1[which(is.infinite(r1[,i])),i] = i0}};rm(i)

jpeg(paste0(pT[2],"density_dbSum.jpeg"), width = 2500, height = 1000, res = 300)
par(mar = c(5,4,2,0)+.1, mfrow = c(1,2))
hist(d.cf[!is.na(d.cf)], freq = T, col = paste0(substr(cBp[1],1,7), "11"), lty = 2, border = cBp[1], main = "", xlab = bquote(italic(d[N]/d[S])~"ratio"), xlim = c(0,2))
hist(d.en[!is.na(d.en)], freq = T, col = paste0(substr(cBp[2],1,7), "11"), lty = 3, border = cBp[2], add = T)

hist(d.cf[!is.na(d.cf)], freq = F, col = paste0(substr(cBp[1],1,7), "11"), lty = 2, border = cBp[1], main = "", xlab = bquote(italic(d[N]/d[S])~"ratio"), xlim = c(0,2))
hist(d.en[!is.na(d.en)], freq = F, col = paste0(substr(cBp[2],1,7), "11"), lty = 3, border = cBp[2], add = T)
legend("topright", legend = c("CF-associated", "Environmental"), fill = cBp[1:2])
invisible(dev.off())

#plot(x = 1:nrow(r0), y = r1[,2]-r1[,5], type = "p", col = cBp[1], pch = 3, cex = .2, xlab = "Sequential order of ORFs", ylab = "log10(dN/dS_CF) - log10(dN/dS_env)") # ylim = c(-1,1)*lIm
#plot(x = 1:nrow(r0), y = r0[,2], type = "l", col = cBp[1], pch = 3, cex = .2, xlab = "Sequential order of ORFs", ylab = "Median dN/dS") # ylim = c(-1,1)*lIm
#lines(x = 1:nrow(r0), y = r0[,5], col = cBp[3], pch = 4, cex = .2)
#legend("topleft", legend = c("CF-associated", "Environmental"), fill = cBp[c(1,3)])
#matplot(r1, type = 'l')
