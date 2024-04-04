#!/bin/bash
# author: ph-u
# script: seqTypingID.r
# desc: concatenate sequences for multi-locus sequence typing (MLST)
# in: source("seqTypingID.r")
# out: NA
# arg: 0
# date: 20240327

# https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=batchSequenceQuery
# https://doi.org/10.12688/wellcomeopenres.14826.1
# https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=query&scheme_id=1
# PA0887-acsA, PA0025-aroE, PA3769-guaA, PA4946-mutL, PA2639-nuoD, PA1770-ppsA, PA0609-trpE
sTIdx = data.frame(PAnum = c("PA0887", "PA0025", "PA3769", "PA4946", "PA2639", "PA1770", "PA0609"), gNam = c("acsA", "aroE", "guaA", "mutL", "nuoD", "ppsA", "trpE"))
sTRef = read.table("../raw/seqTypeRefTable.csv", sep = "\t", header = T)
sTData = read.table("../raw/BIGSdb_2409126_8083931841_34465.csv", sep = "\t", header = T)
sTData = sTData[-grep("[(]",sTData$Locus),]

##### Reallocate sequence type data #####
seqType = as.data.frame(matrix(nr = length(unique(sTData$Contig)), nc = ncol(sTRef)+1))
colnames(seqType) = c("clinical",colnames(sTRef))
seqType$clinical = unique(sTData$Contig)
for(i in 1:nrow(sTData)){
    seqType[which(seqType$clinical==sTData$Contig[i]),which(colnames(seqType)==sTData$Locus[i])] = sTData$Allele[i]
};rm(i)

##### Map sequence type #####
for(i in 1:nrow(seqType)){
    i0 = sTRef
    for(i1 in 2+which(!is.na(seqType[i,-(1:2)]))){
        i0 = i0[which(i0[,i1-1]==seqType[i,i1]),]
        if(nrow(i0)==1){break}
    }
    seqType$ST[i] = paste(unique(i0$ST), collapse = ";")
};rm(i, i0, i1)
