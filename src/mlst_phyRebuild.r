#!/bin/env Rscript
# author: ph-u
# script: mlst_phyRebuild.r
# desc: rebuild and color REALPHY phylogenies
# in: Rscript mlst_phyRebuild.r
# out: res/mlst_REALPHY--*.pdf
# arg: 0
# date: 20250409

source("p_src.r")
outLier = read.tree(paste0(pT[4],"mlst_outliers_polymorphisms_move.phy_phyml_tree.txt"))
outLier$tip.label = read.table(text = gsub("7_gen","7@",gsub("_ASM","@",outLier$tip.label)), sep = "@")[,1]
tipCol = rep("#000000ff",length(outLier$tip.label))
tipCol[-grep("PAO1",outLier$tip.label)] = cBp[as.numeric(as.factor(mEta$sOurce[match(outLier$tip.label[-grep("PAO1",outLier$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))]))]

outLier$edge.length[which(outLier$edge.length==max(outLier$edge.length))] = outLier$edge.length[which(outLier$edge.length==max(outLier$edge.length))]/5

pdf(paste0(pT[2],"mlst_REALPHY--outliers.pdf"), width = 10, height = 20, paper = "a4")
plot(outLier, "u", tip.color = tipCol, cex = .5, main = "Outliers in MLST phylogeny")
legend("topright", legend = levels(as.factor(mEta$sOurce[match(outLier$tip.label[-grep("PAO1",outLier$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))])), fill = cBp[1:(length(unique(tipCol))-1)])
invisible(dev.off())
