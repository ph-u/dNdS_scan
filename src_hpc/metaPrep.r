#!/bin/env Rscript
# author: ph-u
# script: metaPrep.r
# desc: metadata preprocessing
# in: source("metaPrep.r")
# out: NA
# arg: 0
# date: 20240229

mEta = read.csv("../data/metaData.csv", header = T)
mEta[mEta=="missing"]=mEta[mEta==""]=NA

##### f: get rows for necessary category change #####
cHg = function(tErms, tArget, df = mEta){
        i0 = c() 
        for(i in tErms){i0 = unique(c(i0,grep(i,df[,grep("isolation_", colnames(df))])))};rm(i)
        df$sOurce[i0] = tArget
        return(df)
}

mEta$cOuntry = mEta$sOurce = NA
for(i in 1:nrow(mEta)){mEta$cOuntry[i] = ifelse(is.na(mEta[i,grep("geo_loc", colnames(mEta))]),"Unknown",strsplit(mEta[i,grep("geo_loc", colnames(mEta))],":")[[1]][1])};rm(i)
mEta = cHg(c("ystic ", "sputum"), "Cystic fibrosis")
mEta = cHg(c("skin", "toe"), "Skin")
mEta = cHg(c("iver", "water", "soil", "milk", "nvironm", "egg"), "Environmental")
mEta$sOurce[which(is.na(mEta[,grep("isolation_", colnames(mEta))]))] = "Unknown"
mEta$sOurce[is.na(mEta$sOurce)] = "Other infections"

