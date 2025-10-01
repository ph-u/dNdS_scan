#!/bin/env Rscript
# author: ph-u
# script: mlst_stSNP.r
# desc: Get SNP info for each seqType, one gene
# in: Rscript mlst_stSNP.r
# out: data/mlst_stSNP.csv
# arg: 0
# date: 20250411

source("p_src.r")
library(Biostrings);library(msa)
stTarget = read.csv(paste0(pT[2],"mlst_REALPHY--metadata.csv"), header = T)
stList = unique(stTarget$seqType)
oNam = paste0(pT[1],"mlst_stSNP--tmp.fa")

##### set data record format #####
#r0 = as.data.frame(matrix(nrow = nrow(f), ncol = length(stList)))
#row.names(r0) = f$gEne
r1 = vector(mode = "list", length = length(stList))
names(r1) = paste0("ST",stList) # colnames(r0)
for(i in 1:length(r1)){
  r1[[i]] = as.data.frame(matrix(nrow = nrow(f), ncol = choose(sum(stTarget$seqType==stList[i]),2)))
  row.names(r1[[i]]) = f$gEne
};rm(i)

##### get Number of SNP positions for each gene #####
cat(date(),": Start SNP counting\n")
for(i0 in 1:nrow(f)){for(i in 1:length(stList)){
  cat(date(),": Processing gene",i0,"/",nrow(f),"; ST",i,"/",length(stList),"     \r")
  stGene = as.character(read.FASTA(paste0(pT[3],f$fNam[i0]), type = "DNA"))
  s0 = read.table(text = sub("_ASM","@",names(stGene)[-1]), sep = "@")[,1]
  s1 = stTarget[which(stTarget$seqType==stList[i]),1]
  stGene = stGene[which(s0 %in% s1)+1]
  if(length(grep("noHit",names(stGene)))<length(stGene)){
    if(length(grep("noHit",names(stGene)))>0){for(s0 in grep("noHit",names(stGene))){stGene[[s0]] = rep("N",9)}}
    write.FASTA(as.DNAbin(stGene), oNam)
    stGene = msa(readDNAStringSet(oNam), "ClustalW")
    s0 = as.data.frame(strsplit(as.character(stGene),""))
    colnames(s0) = read.table(text = sub("_ASM","@",colnames(s0)), sep = "@")[,1]
    s0 = s0[,s1]
    i3=1;for(i1 in 1:(ncol(s0)-1)){for(i2 in (i1+1):ncol(s0)){r1[[i]][i0,i3] = sum(s0[,i1]!=s0[,i2]);i3 = i3+1}}
#    r0[i0,i] = sum(apply(s0, 1, function(x){length(unique(x))>1}))
}}};rm(i,i0,i1,i2,i3,s0,s1,stGene)
if (file.exists(oNam)) {file.remove(oNam)}

##### export #####
#write.csv(r0[,c(ncol(r0),1:length(stList))], paste0(pT[1],"mlst_stSNP.csv"), row.names = T, quote = F)
for(i in 1:length(stList)){write.csv(r1[[i]],paste0(pT[1],"mlst_stSNP--ST",stList[i],".csv"), row.names = T, quote = F)};rm(i)
