#!/bin/env Rscript
# author: ph-u
# script: reCon.r
# desc: reconstruct residue resolution dN/dS running median
# in: Rscript reCon.r [inFile_db.fa]
# out: data/[strain]_[gene]--reCon.csv
# arg: 1
# date: 20240331 (branch from rDNDS.r)

argv = (commandArgs(T))
#argv = "00_PAO1_107_PA0001_db.fa"
cat(argv,":",date(),"\n")

source("src_dNdS.r"); library(ape)
pT = paste0("../",c("data","binHPC2"),"/")
aRg = list.files(pT[1],argv)
argv = aRg[grep("_db.fa",aRg)]

##### set work env #####
fNam = sub("_db.fa","",argv)
fNam = strsplit(fNam,"--")[[1]]
fNam = c(fNam, sub(paste0("_",fNam[2]),"",fNam[1]))
oNam = paste0(pT[1], paste0(fNam[3:2], collapse = "--"), "--reCon.csv")

rDNDS = read.csv(gsub("reCon","rDNDS",oNam), header = T)

## set rec df
r0.c = c("ntPos","codon","aaRes",paste0("dNdS.",c("min","1Q","2Q","mean","3Q","max")))
f = as.character(read.FASTA(paste0(pT[1], argv[1]), type="DNA"))[[1]]
r0.codon = paste0(f[c(T,F,F)],f[c(F,T,F)],f[c(F,F,T)])
r0 = as.data.frame(matrix(nr=length(r0.codon), nc = length(r0.c)))
colnames(r0) = r0.c
r0$ntPos = seq(1,length(f),3)
r0$codon = r0.codon; rm(r0.codon)
r0$aaRes = c(strsplit(nt2prot(paste0(f, collapse = "")), "")[[1]],"")

##### Calculate running median #####
for(i in 1:nrow(r0)){
    r0.t = rDNDS$dNdS[which(rDNDS$ntStart<=r0$ntPos[i] & rDNDS$ntEnd>=r0$ntPos[i])]
    r0.t[!is.finite(r0.t)] = max(r0.t[is.finite(r0.t)])+1
    r0[i,-(1:3)] = summary(r0.t)
}

##### export #####
write.csv(r0, oNam, row.names = F, quote = F)
cat(argv[1], "residue resolution dN/dS reconstruction done:", date(), "\n")
