#!/bin/env Rscript
# author: ph-u
# script: DataSrc.r
# desc: plot metadata country and sample sources as heatmap
# in: Rscript DataSrc.r
# out: res/DataSrc.pdf
# arg: 0
# date: 20230529

library(lattice)
a = read.csv("../data/metaData.csv", header=T)
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

##### Background distribution #####
pdf("../res/DataSrc.pdf", width = 10, height = 4)
levelplot(table(a[,c("cOuntry","sOurce")]), scales=list(x=list(rot=90)), col.regions = gray(100:0/100), ylab = "Sample source", xlab = "Country")
invisible(dev.off())
