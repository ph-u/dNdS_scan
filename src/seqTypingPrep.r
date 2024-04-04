#!/bin/bash
# author: ph-u
# script: seqTypingPrep.r
# desc: concatenate sequences for multi-locus sequence typing (MLST)
# in: Rscript seqTypingPrep.r
# out: data/seqTypingPrep.fa
# arg: 0
# date: 20240327

# https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=batchSequenceQuery
# https://doi.org/10.12688/wellcomeopenres.14826.1
# https://pubmlst.org/bigsdb?db=pubmlst_paeruginosa_seqdef&page=query&scheme_id=1
# PA0887-acsA, PA0025-aroE, PA3769-guaA, PA4946-mutL, PA2639-nuoD, PA1770-ppsA, PA0609-trpE
source("p_src.r"); library(ape)
x = c("PA0887", "PA0025", "PA3769", "PA4946", "PA2639", "PA1770", "PA0609")
f = list.files(pT[3],"_db.fa")

cat("Collecting sequences for MLST:",date(),"\n")
for(i in 1:length(x)){
    i0 = as.character(read.FASTA(paste0(pT[3],f[grep(paste0(x[i],"_"),f)]), type = "DNA"))
    names(i0) = read.table(text=gsub("_ASM","@",names(i0)), sep = "@")[,1]
    if(i==1){
        r0 = vector(mode = "list", length = length(unique(names(i0))))
        names(r0) = unique(names(i0))
    }; for(i1 in 1:length(r0)){
        cat(i1,"/",length(r0),"(", round(i1/length(r0)*100,2),"% ) in",i,"/",length(x),"(", round(i/length(x)*100,2),"% )",date(),"     \r")
        i2 = grep(names(r0)[i1], names(i0))
        if(length(i2)>0){
            i4=c(0,0);for(i3 in 1:length(i2)){if(length(i0[[i2[i3]]])>i4[2]){i4 = c(i3,length(i0[[i2[i3]]]))}}
            r0[[i1]] = c(r0[[i1]], i0[[i2[i4[1]]]])}
    }
};rm(i, i0, i1, i2, i3, i4)

cat("\nExporting MLST concatenated sequences:",date(),"\n")
write.FASTA(as.DNAbin(r0), paste0(pT[1],"seqTypingPrep.fa"))
cat("Exporting for MLST done:",date(),"\n")
