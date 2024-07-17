#!/bin/env Rscript
# author: ph-u
# script: otherCMD.r
# desc: other commands in data exploration
# in: Rscript: otherCMD.r
# out: NA
# arg: 0
# date: 20240716

source("p_src.r")

##### Number of PAO1-specific genes #####
load(paste0(pT[1],"varTypeBar_circle.rda"))
for(i in 1:length(f.r0)){
    f.r0[[i]]$sRc = names(f.r0)[i]
    f.r0[[i]]$gEne = row.names(f.r0[[i]])
    if(i>1){r0 = rbind(r0, f.r0[[i]])}else{r0 = f.r0[[i]]}
};rm(i)

tHs = 75

cat("\nDetails on the extensive noHit fraction in the PAO1 genome:\n")
r1 = r0[which(r0$noHit>=tHs),c("gEne", "sRc","noHit")]
r2 = as.data.frame(matrix(nr = length(unique(r1$gEne)), nc = length(unique(r1$sRc))))
colnames(r2) = unique(r1$sRc)
row.names(r2) = unique(r1$gEne)
for(i in 1:nrow(r1)){ r2[which(row.names(r2) == r1$gEne[i]),which(colnames(r2) == r1$sRc[i])] = r1$noHit[i] };rm(i)
r2 = round(r2,1)
r2 = r2[which(!is.na(rowSums(r2[,1:2]))),1:2]
#r2[is.na(r2)] = paste0("< ",tHs)
print(r2)

cat("\nDetails on the extensive indel fraction in the PAO1 genome:\n")
r0$INDEL.sel = r0$INDEL.identical + r0$INDEL.SNP
r1 = r0[which(r0$INDEL.sel>=tHs),c("gEne", "sRc","INDEL.sel")]
r3 = as.data.frame(matrix(nr = length(unique(r1$gEne)), nc = length(unique(r1$sRc))))
colnames(r3) = unique(r1$sRc)
row.names(r3) = unique(r1$gEne)
for(i in 1:nrow(r1)){ r3[which(row.names(r3) == r1$gEne[i]),which(colnames(r3) == r1$sRc[i])] = r1$INDEL.sel[i] };rm(i)
r3 = round(r3,1)
r3 = r3[which(!is.na(rowSums(r3[,1:2]))),1:2]
print(r3)

##### Neutral-selected genes #####
cat("\nDetails on neutral-selection along the PAO1 genome, CF vs Environment only:\n")
dNdS.cir = read.csv(paste0(pT[1],"dNdS_median-cir.csv"), header = T)
dC.0 = dNdS.cir[which(dNdS.cir$dNdS==2),-3]
dC.table = as.data.frame.matrix(table(dC.0))
dC.table = dC.table[which(rowSums(dC.table)<ncol(dC.table)),c(1,3)]
dC.table = dC.table[which(rowSums(dC.table)==1),]
cat("\nList of genes neutrally-selected only in CF (",nrow(dC.table[dC.table[,1]==1,]),"):",paste0(paste0(row.names(dC.table[dC.table[,1]==1,]), sep = ", "), collapse = ""),"\n")
cat("\nList of genes neutrally-selected only in the environment (",nrow(dC.table[dC.table[,2]==1,]),"):",paste0(paste0(row.names(dC.table[dC.table[,2]==1,]), sep = ", "), collapse = ""),"\n")
