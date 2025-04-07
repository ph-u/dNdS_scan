#!/bin/env Rscript
# author: ph-u
# script: pathogenwatch_seqType.r
# desc: Plot sequence type distribution
# in: Rscript pathogenwatch_seqType.r
# out: ../res/pathogenwatch_seqType.pdf
# arg: 0
# date: 20250407

cf = read.csv("../raw/pathogenwatch--8qzi0c5dot7p-paipcdinitial_cfnenv-typing.csv", header = T)
aLl = read.csv("../raw/pathogenwatch--j8j8npqh4cd4-paipcdinitial_all-typing.csv", header = T)

aLl$cfenv = cf$MLST.ST..Pseudomonas.aeruginosa..PubMLST..[match(aLl$NAME, cf$NAME)]
a0 = aLl[which(!is.na(aLl$cfenv)),]
all(a0$MLST.ST..Pseudomonas.aeruginosa..PubMLST..==a0$cfenv) # TRUE

##### plot seq-type distribution table #####
pdf("../res/pathogenwatch_seqType.pdf", width = 10, height = 20, paper = "a4")
par(mfrow = c(2,1))
x = table(cf$MLST.ST..Pseudomonas.aeruginosa..PubMLST..)
plot(x, ylab = "sequence type frequency", xlab = "sequence type", main = paste0("CF/environmental isolates: n = ",sum(x)), xaxt = "n", ylim = c(0,ceiling(max(x)/10)*10))
x = table(aLl$MLST.ST..Pseudomonas.aeruginosa..PubMLST..)
plot(x, ylab = "sequence type frequency", xlab = "sequence type", main = paste0("all isolates: n = ",sum(x)), xaxt = "n", ylim = c(0,ceiling(max(x)/10)*10))
invisible(dev.off())
