#!/bin/env Rscript
# author: ph-u
# script: dNdS_convert.r
# desc: extract source-specific dNdS mapping information
# in: Rscript dNdS_convert.r [PAnum]
# out: data/[strain]_[PAnum]_[env]--reCon.csv
# arg: 1
# date: 20240405 (mod from binHPC/reCon.r)

argv = (commandArgs(T))
#argv = "PA5221"
pRe = "PAO1_107_"

source("p_src.r"); library(ape)
rDNDS = read.csv(paste0(pT[3], pRe, argv[1], "--rDNDS.csv"), header = T)
r0.c = c("ntPos","codon","aaRes",paste0("dNdS.",c("min","1Q","2Q","mean","3Q","max")))
f = as.character(read.FASTA(paste0(pT[3], "00_", pRe, argv[1], "_db.fa"), type="DNA"))[[1]]
r0.codon = paste0(f[c(T,F,F)],f[c(F,T,F)],f[c(F,F,T)])

sRc = unique(mEta$sOurce)
for(i in 1:length(sRc)){
    r1 = rDNDS[which(read.table(text = gsub("_ASM","@",rDNDS$clinical), sep = "@")[,1] %in% mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==sRc[i])]),]

    if(nrow(r1)>1){
        r0 = as.data.frame(matrix(nr=length(r0.codon), nc = length(r0.c)))
        colnames(r0) = r0.c
        r0$ntPos = seq(1,length(f),3)
        r0$codon = r0.codon
        r0$aaRes = c(strsplit(nt2prot(paste0(f, collapse = "")), "")[[1]],"")

##### Calculate running median #####
        for(i0 in 1:nrow(r0)){
            r0.t = r1$dNdS[which(r1$ntStart<=r0$ntPos[i0] & r1$ntEnd>=r0$ntPos[i0])]
            r0.t[!is.finite(r0.t)] = max(r0.t[is.finite(r0.t)])+1
            r0[i0,-(1:3)] = summary(r0.t)
        }


        write.csv(r0, paste0(pT[1], pRe, argv[1], "_", gsub(" ", "-", sRc[i]), "--reCon.csv"), row.names = F, quote = F)
}};rm(i)
