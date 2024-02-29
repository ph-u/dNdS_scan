#!/bin/env Rscript
# author: ph-u
# script: varTypeMap.r
# desc: plot dN/dS variation types as heatmap
# in: Rscript varTypeMap.r [srcFile] [strain] [gene]
# out: res/04_[strain]_[gene]_vtMap{Src,Loc}.pdf
# arg: 3
# date: 20240229

#argv = (commandArgs(T))
argv = c("00_PAO1_107_PA2668.csv", "PAO1_107", "PA2668")
pT = paste0("../",c("data","res"),"/")
library(lattice)

dNdS = read.csv(paste0(pT[1],argv[1]), header = T)
source("metaPrep.r")

##### Sequence variation type extraction #####
dNdS$varType = dNdS$sampleEnv = dNdS$cOuntry = NA
for(i in grep("%", dNdS$locus)){
    dNdS$varType[i] = strsplit(dNdS$locus[i],"%")[[1]][1]
};rm(i)
dNdS$varType[is.na(dNdS$varType)] = "SNP"
for(i in 1:nrow(mEta)){
    i0 = grep(mEta$assemblyInfo.genbankAssmAccession[i], dNdS$clinical)
    dNdS$sampleEnv[i0] = mEta$sOurce[i]
    dNdS$cOuntry[i0] = mEta$cOuntry[i]
};rm(i, i0)

##### plot #####
for(i in 1:2){
    if(i<2){
        i0 = c("Src", "sampleEnv", "Sample Source")
    }else{
        i0 = c("Loc", "cOuntry", "Sampling Country")
    }
    pdf(paste0(pT[2],"04_",argv[2],"_",argv[3],"vtMap",i0[1],".pdf"), width = 10, height = 4)
    print(levelplot(table(dNdS[,c(i0[2], "varType")]), scales=list(x=list(rot=90)), col.regions = gray(100:0/100), ylab = "Seq Variation Type", xlab = i0[3]))
    invisible(dev.off())
};rm(i)
