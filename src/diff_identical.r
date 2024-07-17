#!/bin/env Rscript
# author: ph-u
# script: diff_identical.r
# desc: Get gene list which dN/dS of CF different from non-CF sources
# in: Rscript diff_identical.r
# out: NA
# arg: 0
# date: 20240422

source("p_src.r")
f = f[match(protList, f[,2]),] # filter only proteomics-confirmed fraction

##### Calculate % identical for each gene among sampling sources #####
r0 = as.data.frame(matrix(0, nr = nrow(f), nc = length(unique(mEta$sOurce))))
row.names(r0) = f[,2]; colnames(r0) = unique(mEta$sOurce)

cat("Summarizing Identical gene %:",date(),"\n")
for(i in 1:nrow(f)){
    cat("Processing gene",f[i,2],":",date(),"     \r")
    f.0 = read.csv(paste0(pT[3],"PAO1_107_",f[i,2],"--dbSum.csv"), header = T)
    f.0$sOurce = mEta$sOurce[match(mEta$assemblyInfo.genbankAssmAccession, read.table(text = gsub("_ASM","@",f.0$clinical), sep = "@")[,1])]
    for(i0 in grep(":",f.0$varType)){
        i1 = strsplit(f.0$varType[i0], ":")[[1]]
        f.0$varType[i0] = paste(i1[-c(3:6)], collapse = ":")}
    if(any(f.0$varType=="identical")){
        r.0 = table(f.0[,c("varType","sOurce")])
        r.0 = r.0["identical",]/colSums(r.0)
        r0[i,] = r.0[match(colnames(r0), names(r.0))]}
};rm(i,i0,i1,f.0, r.0);cat("\nSummary done:",date(),"\n")
#r0 = r0[rowSums(r0)>0,]*100

##### Pairwise Analysis #####
r0.L = data.frame(gene=row.names(r0), src=rep(colnames(r0), each = nrow(r0)), perc=unname(unlist(r0)))
print(pairwise.wilcox.test(r0.L$perc, r0.L$src, p.adjust = "bonf"))

cOlr0 = data.frame(diff = apply(r0, 1, function(x) diff(range(x))), cOl="#00000005")
cOlr0[cOlr0.0 <- rev(order(cOlr0$diff))[1:20],2] = "#000000ff"
cOlr0[which(row.names(r0)=="PA2185"),2] = "#ff00ffff"
r0 = r0[,order(colnames(r0))]
r0.0 = r0[cOlr0.0,]

print(pairwise.wilcox.test(r0.L$perc[which(r0.L$gene %in% row.names(r0.0))], r0.L$src[which(r0.L$gene %in% row.names(r0.0))], p.adjust = "bonf"))

pdf(paste0(pT[2],"idPerc_level.pdf"), width=6, height=7)
par(mar = c(3,4,0,1)+.1)
matplot(t(r0)*100, type="b", cex=.1, pch=18, lty = 1, col = cOlr0$cOl, xaxt = "n", ylab = "Clinical Isolates with identical sequence with PAO1 (%)")
axis(1, at = 1:ncol(r0), labels = gsub(" ","\n",colnames(r0)), padj = .3)
text(x=rep(1:ncol(r0.0), each = nrow(r0.0)), y=unlist(r0.0)*100, cex = .5, labels=row.names(r0.0))
invisible(dev.off())
