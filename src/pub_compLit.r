#!/bin/env Rscript
# author: ph-u
# script: pub_compLit.r
# desc: compare result with literature
# in: Rscript pub_compLit.r
# out: res/p_*--compLit.{pdf,cvs}
# arg: 0
# date: 20250327

source("p_src.r")
d0 = read.csv(paste0(pT[1],"p_conservationTables.csv"), header = T)
d1 = read.table(paste0(pT[4],"pnas.1419677112.sd04.csv"), skip = 1, sep = "\t", comment.char = "", quote = "", header = T)
f.kmeansDF = read.csv(paste0(pT[2], "string_diff.csv"), header = T)

d1.m = match(d0$gene, d1$PAO1.Locus.ID)
# d1 = d1[which(d1$PAO1.Locus.ID %in% d0$gene),] # 5569 ORFs
# d0[which(!(d0$gene %in% d1$PAO1.Locus.ID)),] # 17 ORFs not shared between d0 & d1
d0$PAO1.log2SpuExp = d1$Log2.Fold.Change.Mutant.Abundance..Sputum.vs..Expected..1[d1.m]
d0$PA14.log2SpuExp = d1$Log2.Fold.Change.Mutant.Abundance..Sputum.vs..Expected.[d1.m]
d0$volCol = ifelse(d0$volcano.cOlor=="Grey", "#00000011", ifelse(d0$volcano.cOlor=="Blue", "#00001133", ifelse(d0$volcano.cOlor=="Red", "#FF0000FF", "#00FF00FF")))
d0.b = d0
d0.m = range(d0$log2FC[is.finite(d0$log2FC)], na.rm = T)
d0$log2FC[d0$log2FC==-Inf] = d0.m[1]-1
d0$log2FC[d0$log2FC==Inf] = d0.m[2]+1

##### O: plot log2FC comparison graph #####
pdf(paste0(pT[2],"p_log2FC.1--compLit.pdf"), width = 10, height = 20, paper = "a4")
#plot(y = d1$Log2.Fold.Change.Mutant.Abundance..Sputum.vs..Expected., x = d1$Log2.Fold.Change.Mutant.Abundance..Sputum.vs..Expected..1, pch = 20, xlab = "PAO1 log2(sputum vs expected)", ylab = "PA14 log2(sputum vs expected)")
plot(y = d0$PA14.log2SpuExp, x = d0$PAO1.log2SpuExp, pch = 20, xlab = "PAO1 log2(sputum vs expected)", ylab = "PA14 log2(sputum vs expected)", main = "x = Literature-PAO1, y = Literature-PA14", col = d0$volCol)
abline(h = 0, v = 0)

#plot(x = d0$log2FC, y = d1$Log2.Fold.Change.Mutant.Abundance..Sputum.vs..Expected..1, pch = 20, ylab = "PAO1 log2(sputum vs expected)", xlab = "Selection pressure log2(CF vs env)")
plot(x = d0$log2FC, y = d0$PAO1.log2SpuExp, pch = 20, ylab = "PAO1 log2(sputum vs expected)", xlab = "Selection pressure log2(CF vs env)", main = "x = dN/dS, y = Literature-PAO1", col = d0$volCol)
abline(h = c(0,-6.5), v = 0)

#plot(y = d1$Log2.Fold.Change.Mutant.Abundance..Sputum.vs..Expected., x = d0$log2FC, pch = 20, xlab = "Selection pressure log2(CF vs env)", ylab = "PA14 log2(sputum vs expected)")
plot(y = d0$PA14.log2SpuExp, x = d0$log2FC, pch = 20, xlab = "Selection pressure log2(CF vs env)", ylab = "PA14 log2(sputum vs expected)", main = "x = dN/dS, y = Literature-PA14", col = d0$volCol)
abline(h = c(0,-6.5), v = 0)
invisible(dev.off())

##### O: PCA #####
d00.pca = pcaSwap(d0[which(!is.na(d0[,13]) & !is.na(d0[,14])),c(7,13,14)])
d9 = pca(d0[which(!is.na(d0[,13]) & !is.na(d0[,14])),c(7,13,14)], method = "ppca")
d00.plot = ggbiplot(d00.pca)
#ggsave(paste0(pT[2],"p_log2FC.2--compLit.pdf"), plot = pcaLAB(d00.plot,round(d9@R2*100,2)))

##### extract high log2FC shown in dN/dS vs Lit-PAO1/PA14 plot #####
d0 = d0.b
d2.PAO1 = d0[which((d0$volcano.cOlor %in% c("Red", "Green")) & (d0$PAO1.log2SpuExp < -6.5)),] # nrow = 46
d2.PA14 = d0[which((d0$volcano.cOlor %in% c("Red", "Green")) & (d0$PA14.log2SpuExp < -6.5)),] # nrow = 24
write.csv(d2.PAO1[,c(1,2,7,10,11,13,14,12)], paste0(pT[2],"p_litPAO1--compLit.csv"), quote = F, row.names = F)
write.csv(d2.PA14[,c(1,2,7,10,11,13,14,12)], paste0(pT[2],"p_litPA14--compLit.csv"), quote = F, row.names = F)

##### match STRING consernsus clusters with volcano plot extremes #####
d3 = d0[which((d0$log2FC==Inf) & (d0$log10P>-log10(.1))),] # 111 volcano right-most ORFs (has 12 Green?)
d4 = d0[which((d0$log2FC==-Inf) & (d0$log10P>-log10(.1))),] # 19 volcano left-most ORFs (has 4 Green?)

d0$STRING.cluster = f.kmeansDF$cID[match(d0$gene,f.kmeansDF$PAnum)]

d3 = d0[which((d0$log2FC==Inf) & (d0$volcano.cOlor=="Red")),] # 99 volcano right-most ORFs
d3 = d3[,c(1,2,10,11,16)]
d3 = d3[rev(order(d3$log10P)),]

d4 = d0[which((d0$log2FC==-Inf) & (d0$volcano.cOlor=="Red")),] # 15 volcano left-most ORFs
d4 = d4[,c(1,2,10,11,16)]
d4 = d4[rev(order(d4$log10P)),]

d3[is.na(d3)] = ""
d4[is.na(d4)] = ""
write.csv(d3, paste0(pT[2],"p_volcanoRightEx--compLit.csv"), row.names = F, quote = F)
write.csv(d4, paste0(pT[2],"p_volcanoLeftEx--compLit.csv"), row.names = F, quote = F)
d3[d3==""] = NA
d4[d4==""] = NA

## Cluster enrichment evaluation
##### f: enrichment score #####
enRich = function(dfTest, dfRef = as.data.frame(table(f.kmeansDF$cID))){
  sCore = rep(0, nrow(dfRef))
  names(sCore) = dfRef$Var1
  xTab = as.data.frame(table(dfTest$STRING.cluster))
  xTab[,1] = as.character(xTab[,1])
  for(i in 1:length(sCore)){ if(i %in% xTab$Var1){
    sCore[i] = xTab$Freq[xTab$Var1==names(sCore)[i]]
  }}
  return(sCore / (dfRef$Freq*nrow(dfTest)/sum(dfRef$Freq)))} # Observed vs Expected

enRich(d3)
enRich(d4)
