#!/bin/env Rscript
# author: ph-u
# script: PAO1_geneGC.r
# desc: calculate gene-wise GC content
# in: Rscript PAO1_geneGC.r
# out: res/PAO1_geneGC.pdf
# arg: 0
# date: 20240823

source("p_src.r")
x = as.character(read.FASTA(paste0(pT[4],"PAO1_107.fa"), type="DNA"))
f$GC = NA
for(i in 1:nrow(f)){f$GC[i] = GC.content(as.DNAbin(x[[i]]))};rm(i)

xInt = grep("PA2125|PA2384", f$gEne)

pdf(paste0(pT[2],"PAO1_geneGC.pdf"))
plot(x=1:nrow(f),y=f$GC*100, type = "l", xlab = "Gene position", ylab = "GC content (%)", col = cBp[5])
abline(h = c(GC.content(as.DNAbin(x)),median(f$GC))*100, v = c(1900,2600), lty = 3)
abline(h = GC.content(as.DNAbin(x[xInt[1]:xInt[2]]))*100, v = xInt, lty = 2, col = cBp[2])
invisible(dev.off())

print(wilcox.test(f$GC[xInt[1]:xInt[2]], f$GC[1900:2600]))
