#!/bin/env Rscript
# author: ph-u
# script: cf_aaAlign.r
# desc: extract protein sequences for CF only strains
# in: Rscript cf_aaAlign.r
# out: data/*-cfOnly.fa
# arg: 0
# date: 20240210

library(ape);library(msa)
pT = "../data/"
a = read.csv(paste0(pT,"metaData.csv"), header=T)
a[a=="missing"]=a[a==""]=NA

##### f: get rows for necessary category change #####
cHg = function(tErms, tArget, df = a){
	i0 = c()
	for(i in tErms){i0 = unique(c(i0,grep(i,df[,grep("isolation_", colnames(df))])))};rm(i)
	df$sOurce[i0] = tArget
	return(df)}

a$cOuntry = a$sOurce = NA
for(i in 1:nrow(a)){a$cOuntry[i] = ifelse(is.na(a[i,grep("geo_loc", colnames(a))]),"Unknown",strsplit(a[i,grep("geo_loc", colnames(a))],":")[[1]][1])};rm(i)
a = cHg(c("ystic ", "sputum"), "Cystic fibrosis")
a = cHg(c("skin", "toe"), "Skin")
a = cHg(c("iver", "water", "soil", "milk", "nvironm", "egg"), "Environmental")
a$sOurce[which(is.na(a[,grep("isolation_", colnames(a))]))] = "Unknown"
a$sOurce[is.na(a$sOurce)] = "Other infections"

x = a$assemblyInfo.genbankAssmAccession[which(a$sOurce=="Cystic fibrosis")]
for(i in list.files("../data","aaAlign.fa")){
    x0 = vector(mode="list", length = length(x))
    a0 = as.character(read.FASTA(paste0(pT,i), type = "AA"))
    for(i0 in 1:length(x)){
        x1 = grep(x[i0], names(a0))
        x0[[i0]] = a0[[x1]]
        names(x0)[i0] = names(a0)[x1]
    }
    write.FASTA(as.AAbin(x0), paste0(pT,gsub("aaAlign.fa","aaAlign-cfOnly.fa",i)))
    a0 = msa(readAAStringSet(paste0(pT,gsub("aaAlign.fa","aaAlign-cfOnly.fa",i))), "ClustalOmega")
    writeXStringSet(unmasked(a0), file=paste0(pT,gsub("aaAlign.fa","aaAlign-cfOnly.fa",i)))
};rm(i,i0,x0,a0,x1)
