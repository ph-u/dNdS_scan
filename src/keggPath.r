#!/bin/env Rscript
# author: ph-u
# script: keggPath.r
# desc: statistical test on habitat-specific selection pressure
# in: Rscript keggPath.r
# out: NA
# arg: 0
# date: 20250610

##### metabolic pathway #####
a = read.table("../raw/keggPath.csv", sep = "\t", header = T)
a[a=="-"] = 0
for(i in 3:ncol(a)){a[,i] = as.numeric(a[,i])};rm(i)
#w0 = wilcox.test(c(a$CF.Within.pathway,a$env.Within.pathway),c(a$CF.Associated,a$env.Associated))

a0.cf = a[,-c(5,6)];a0.cf$src = "CF-associated"
a0.ev = a[,-c(3,4)];a0.ev$src = "Environmental"
colnames(a0.cf)[3:4] = colnames(a0.ev)[3:4] = c("Within.pathway","Associated")
a0 = rbind(a0.cf,a0.ev);rm(a0.cf,a0.ev)
w0 = wilcox.test(a0$Within.pathway,a0$Associated)

a.bx = data.frame(src=rep(a0$src,2), path=rep(c("Within.pathway","Associated"), each = nrow(a0)), count = c(a0[,3],a0[,4]))

#pdf("../res/keggPath.pdf")
#par(mar = c(4,4,0,0)+.1, cex = 2)
#boxplot(count ~ path, data = a.bx, col = "#00000000", ylab = "Count", xlab = "ORF-metabolism relationship", pch = 4)
#segments(x0 = 1, x1 = 2, y0 = 6, y1 = 6)
#text(x = 1.5, y = 7.5, label = paste0("Wilcox test\nW = ",w0$statistic, ", p = ",round(w0$p.value,4)))
#invisible(dev.off())

##### metabolites #####
m0 = read.table("../raw/kegg_metabolites.csv", sep = "\t", header = T)
colnames(m0)[(-1:0)+ncol(m0)] = c("CF", "Environmental")
w0 = wilcox.test(m0$CF,m0$Environmental)
