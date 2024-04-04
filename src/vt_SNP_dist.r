#!/bin/env Rscript
# author: ph-u
# script: vt_SNP_dist.r
# desc: extract SNP extremes
# in: Rscript vt_SNP_dist.r
# out: NA
# arg: 0
# date: 20240321

source("p_src.r")
teS = read.csv(paste0(pT[2],"vt_Bar.csv"), row.names = 1, header = T)
teS$SNP = teS$SNP + teS$INDEL.SNP
snp.Lst = row.names(teS)[rev(order(teS$SNP[which(teS$SNP >= 99)]))] # top 1% SNP-dominated genes
#f = f[which(f$gEne %in% row.names(teS)[teS$SNP >= 50]),]

rEs.c = c("conserved", "neutral", "variable")
rEs = as.data.frame(matrix(nr = length(snp.Lst), nc = length(rEs.c)))
colnames(rEs) = rEs.c; rm(rEs.c)
row.names(rEs) = snp.Lst
cat("Scan dN/dS distribuion for sequence segments:",date(),"\n")
for(i in 1:nrow(rEs)){
    cat(i,"/",nrow(rEs),":",date(), "     \r")
    i0 = read.csv(paste0(pT[3],"PAO1_107_",row.names(rEs)[i],"--dbSum.csv"), header = T)
    i0 = i0[grep("SNP", i0$varType),]
    if(sum(is.finite(i0$dNdS))>0){
        i0$dNdS[!is.finite(i0$dNdS)] = max(i0$dNdS[is.finite(i0$dNdS)])+1
        rEs[i,] = c(sum(i0$dNdS <= .5), sum(abs(i0$dNdS-1) <= .05), sum(i0$dNdS >= 2))/nrow(i0)*100
    }else{rEs[i,] = rep(NA,3)}
};rm(i, i0)

cat("\nWriting result:",date(),"\n")
write.csv(rEs, paste0(pT[1],"vt_SNP_dist.csv"), row.names = T, quote = F)
cat("Exported:",date(),"\n")

nLst = 20
topList = as.data.frame(matrix(nr = nLst, nc = ncol(rEs)))
colnames(topList) = colnames(rEs)
for(i in 1:ncol(rEs)){
    i0 = row.names(rEs)[rev(order(rEs[!is.na(rEs[,i]),i]))[1:nLst]]
    topList[,i] = i0[order(i0)]
    cat("Top", nLst, "genes with", colnames(rEs)[i], "selection regions:", paste(i0, collapse = ", "), "\n\n")
};rm(i, i0)
write.csv(topList, paste0(pT[2],"vt_SNP_topList.csv"), row.names = F, quote = F)
