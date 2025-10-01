#!/bin/env Rscript
# author: ph-u
# script: pathogenwatch_seqType.r
# desc: Plot sequence type distribution
# in: Rscript pathogenwatch_seqType.r
# out: ../res/pathogenwatch_seqType.pdf
# arg: 0
# date: 20250407

source("p_src.r")
cf = read.csv(paste0(pT[4],"pathogenwatch--8qzi0c5dot7p-paipcdinitial_cfnenv-typing.csv"), header = T)
cf$sOurce = mEta$sOurce[match(cf$NAME,mEta$assemblyInfo.genbankAssmAccession)]
cf$cOuntry = mEta$cOuntry[match(cf$NAME,mEta$assemblyInfo.genbankAssmAccession)]
aLl = read.csv(paste0(pT[4],"pathogenwatch--j8j8npqh4cd4-paipcdinitial_all-typing.csv"), header = T)
aLl$sOurce = mEta$sOurce[match(aLl$NAME,mEta$assemblyInfo.genbankAssmAccession)]
aLl$cOuntry = mEta$cOuntry[match(aLl$NAME,mEta$assemblyInfo.genbankAssmAccession)]

aLl$cfenv = cf$MLST.ST..Pseudomonas.aeruginosa..PubMLST..[match(aLl$NAME, cf$NAME)]
a0 = aLl[which(!is.na(aLl$cfenv)),]
all(a0$MLST.ST..Pseudomonas.aeruginosa..PubMLST..==a0$cfenv) # TRUE

##### plot seq-type distribution table #####
pdf(paste0(pT[2],"pathogenwatch_seqType.pdf"), width = 10, height = 20, paper = "a4")
par(mfrow = c(2,1))
x = table(cf[,c(2,5)])
barplot(t(x), ylab = "sequence type frequency", xlab = "sequence type", main = paste0("CF/environmental isolates: n = ",sum(x)), xaxt = "n", ylim = c(0,ceiling(max(x)/10)*10), col=cBp[1:ncol(x)], border = "#00000000")
text(x = (1:nrow(x))*1.2, y = rowSums(x)+1.2, labels = ifelse(rowSums(x)>10,ifelse(nchar(rownames(x))>5,"",rownames(x)),""), srt = 90)
legend("topright", legend = colnames(x), fill = cBp[1:ncol(x)])

x = cf[which(cf$MLST.ST..Pseudomonas.aeruginosa..PubMLST.. %in% rownames(x)[which(rowSums(x)>10)]),-c(3:4)]
for(i in ncol(x):2){x = x[order(x[,i]),]};rm(i)
x$originalSampleType = mEta$assemblyInfo.biosample.attributes.sample_type[match(x$NAME,mEta$assemblyInfo.genbankAssmAccession)]
x$originalIsolationSource = gsub(",",";",mEta$assemblyInfo.biosample.attributes.isolation_source[match(x$NAME,mEta$assemblyInfo.genbankAssmAccession)])
write.csv(x, paste0(pT[2],"pathogenwatch_seqType.csv"), quote = F, row.names = F)

x = table(aLl[,c(2,5)])
barplot(t(x), ylab = "sequence type frequency", xlab = "sequence type", main = paste0("all isolates: n = ",sum(x)), xaxt = "n", ylim = c(0,ceiling(max(rowSums(x))/10)*10), col=cBp[1:ncol(x)], border = "#00000000")
text(x = (1:nrow(x))*1.2, y = rowSums(x)+2, labels = ifelse(rowSums(x)>10,ifelse(nchar(rownames(x))>5,"",rownames(x)),""), srt = 90)
legend("topright", legend = colnames(x), fill = cBp[1:ncol(x)])
invisible(dev.off())

jpeg(paste0(pT[2],"pathogenwatch_seqType.jpeg"), width = 1000, height = 700, res = 100)
x = table(cf[,c(2,5)])
barplot(t(x), ylab = "sequence type frequency", xlab = "sequence type", main = paste0("CF/environmental isolates: n = ",sum(x)), xaxt = "n", ylim = c(0,ceiling(max(x)/10)*10), col=cBp[(1:ncol(x))+2], border = "#00000000")
text(x = (1:nrow(x))*1.2, y = rowSums(x)+1.2, labels = ifelse(rowSums(x)>10,ifelse(nchar(rownames(x))>5,"",paste0("  ST",rownames(x))),""), srt = 90)
legend("topright", legend = colnames(x), fill = cBp[(1:ncol(x))+2])
invisible(dev.off())
