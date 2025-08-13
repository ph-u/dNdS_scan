#!/bin/env Rscript
# author: ph-u
# script: mlst_phyRebuild.r
# desc: rebuild and color REALPHY phylogenies
# in: Rscript mlst_phyRebuild.r
# out: res/mlst_REALPHY--*.pdf
# arg: 0
# date: 20250409

source("p_src.r")
oCex = c(1500,3500,.0056)
outLier = read.tree(paste0(pT[4],"mlst_outliers_polymorphisms_move.phy_phyml_tree.txt"))
outLier$tip.label = read.table(text = gsub("7_gen","7@",gsub("_ASM","@",outLier$tip.label)), sep = "@")[,1]
tipCol = tipCol0 = rep("#000000ff",length(outLier$tip.label))
tipCol[-grep("PAO1",outLier$tip.label)] = cBp[as.numeric(as.factor(mEta$sOurce[match(outLier$tip.label[-grep("PAO1",outLier$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))]))]
tipCol0[-grep("PAO1",outLier$tip.label)] = cBp[as.numeric(as.factor(mEta$cOuntry[match(outLier$tip.label[-grep("PAO1",outLier$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))]))]

#outLier$edge.length[which(outLier$edge.length==max(outLier$edge.length))] = outLier$edge.length[which(outLier$edge.length==max(outLier$edge.length))]/5

jpeg(paste0(pT[2],"mlst_REALPHY--outliers_src.jpeg"), width = oCex[1], height = oCex[1], res = 300)
#pdf(paste0(pT[2],"mlst_REALPHY--outliers.pdf"), width = 10, height = 20, paper = "a4")
plot(outLier, "u", tip.color = tipCol, cex = .5, main = "Outliers in MLST phylogeny")
legend("topright", legend = levels(as.factor(mEta$sOurce[match(outLier$tip.label[-grep("PAO1",outLier$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))])), fill = cBp[1:(length(unique(tipCol))-1)])
invisible(dev.off())

jpeg(paste0(pT[2],"mlst_REALPHY--outliers_geo.jpeg"), width = oCex[1], height = oCex[1], res = 300)
plot(outLier, "u", tip.color = tipCol0, cex = .5, main = "Outliers in MLST phylogeny")
legend("topright", legend = levels(as.factor(mEta$cOuntry[match(outLier$tip.label[-grep("PAO1",outLier$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))])), fill = cBp[1:(length(unique(tipCol0))-1)])
invisible(dev.off())


stDist = read.tree(paste0(pT[4],"realphy_seqTypeDistri_polymorphisms_move.phy_phyml_tree.txt"))
stData = read.csv(paste0(pT[4],"pathogenwatch--8qzi0c5dot7p-paipcdinitial_cfnenv-typing.csv"), header = T)
stDist$tip.label = read.table(text = gsub("7_gen","7@",gsub("_ASM","@",stDist$tip.label)), sep = "@")[,1]
tipMtx = matrix(1,nrow = length(stDist$tip.label), ncol = 3)
row.names(tipMtx) = stDist$tip.label
tipCol = tipCol0 = tipCol1 = rep("#00000000",length(stDist$tip.label))
tipCol[-grep("PAO1",stDist$tip.label)] = cBp[as.numeric(as.factor(mEta$sOurce[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))]))]
tipCol0[-grep("PAO1",stDist$tip.label)] = cBp[as.numeric(as.factor(mEta$cOuntry[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))]))]
tipCol1[-grep("PAO1",stDist$tip.label)] = cBp[as.numeric(as.factor(mEta$seqType[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))]))]
#stDist$tip.label[-grep("PAO1",stDist$tip.label)] = paste0("ST",stData[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",stData$NAME)),2])
#stDist$tip.label[-grep("PAO1",stDist$tip.label)] = paste0("ST",stData[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",stData$NAME)),2]," ",stDist$tip.label[-grep("PAO1",stDist$tip.label)]," ST",stData[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",stData$NAME)),2])
tCol = list(tipCol0,tipCol,tipCol1)
stDist = root(stDist,outgroup = grep("PAO1",stDist$tip.label))

jpeg(paste0(pT[2],"mlst_REALPHY--stDist.jpeg"), width = oCex[2], height = oCex[2], res = 300)
par(mai = rep(1,4)+.1, xpd = T)
#plot(stDist, tip.color = tipCol0, show.tip.label = T, cex = .5, no.margin = T)
#for(i in 1:ncol(tipMtx)){phydataplot(tipMtx[,i]/2000, stDist, style = "bars", border = "#00000000", offset = i/2000-.005, col = tCol[[i]])}
plot(stDist, "f", x.lim = c(-oCex[3],oCex[3]), y.lim = c(-oCex[3],oCex[3]), tip.color = "#000000ff", show.tip.label = T, cex = .5, no.margin = T)
for(i in 1:ncol(tipMtx)){ring(tipMtx[,i]/2000, stDist, style = "ring", offset = i/2000+3e-4, col = tCol[[i]])}
invisible(dev.off())

#jpeg(paste0(pT[2],"mlst_REALPHY--stDist_src.jpeg"), width = oCex[2], height = oCex[2], res = 300)
#plot(stDist, "f", tip.color = tipCol, cex = .5, main = "SeqType in MLST phylogeny")
#invisible(dev.off())

#jpeg(paste0(pT[2],"mlst_REALPHY--stDist_geo.jpeg"), width = oCex[2], height = oCex[2], res = 300)
#plot(stDist, "f", tip.color = tipCol0, cex = .5, main = "SeqType in MLST phylogeny")
#invisible(dev.off())

jpeg(paste0(pT[2],"mlst_REALPHY--stDist_L.jpeg"), width = oCex, height = oCex, res = 300)
par(mar = rep(0,4)+.1)
plot(0,0, col = "#00000000", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("topright", legend = levels(as.factor(mEta$sOurce[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))])), fill = cBp[1:(length(unique(tipCol))-1)])
legend("topleft", legend = levels(as.factor(mEta$cOuntry[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))])), fill = cBp[1:(length(unique(tipCol0))-1)])
legend("bottomright", legend = paste0("ST",levels(as.factor(mEta$seqType[match(stDist$tip.label[-grep("PAO1",stDist$tip.label)],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))]))), fill = cBp[1:(length(unique(tipCol1))-1)])
#legend("topleft", legend = levels(as.factor(mEta$sOurce[match(read.table(text = stDist$tip.label[-grep("PAO1",stDist$tip.label)], sep = " ")[,2],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))])), fill = cBp[1:(length(unique(tipCol))-1)])
#legend("topright", legend = levels(as.factor(mEta$cOuntry[match(read.table(text = stDist$tip.label[-grep("PAO1",stDist$tip.label)], sep = " ")[,2],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))])), fill = cBp[1:(length(unique(tipCol0))-1)])
#legend("bottomleft", legend = paste0("ST",levels(as.factor(mEta$seqType[match(read.table(text = stDist$tip.label[-grep("PAO1",stDist$tip.label)], sep = " ")[,2],gsub("[.]","_",mEta$assemblyInfo.genbankAssmAccession))]))), fill = cBp[1:(length(unique(tipCol1))-1)])
invisible(dev.off())

##### Extract metadata into from same seqeunce type, shared between CF & env within same country #####
m0 = mEta[which(mEta$assemblyInfo.genbankAssmAccession %in% paste0("GCA_00",c(437542,383338,437416,437449,437448,437450,437552,397633,397637,383988,383996,383896,383887,383889,383881),"5.1")),c(3,16:18,13,15)]
m0[,ncol(m0)] = gsub(",",";",m0[,ncol(m0)])
write.csv(m0[order(m0$seqType),], paste0(pT[2],"mlst_REALPHY--metadata.csv"), quote = F, row.names = F)
