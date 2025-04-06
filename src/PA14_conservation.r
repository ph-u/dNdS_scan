#!/bin/env Rscript
# author: ph-u
# script: PA14_conservation.r
# desc: conservation categorization for PA14 orthologs of the PAO1 wrong sequences
# in: Rscript PA14_conservation.r
# out: stdout
# arg: 0
# date: 20230323

##### get tHs0 from dNdS_median-cir.r #####
source("p_src.r")
source("p_metabolism_PAO1.r")
dNdS.med = read.csv(paste0(pT[1],"dNdS_median-PCA.csv"), header = T, row.names = 1)[,c(2,4,3,5)]
colnames(dNdS.med) = paste(1:ncol(dNdS.med),colnames(dNdS.med), sep = ".")
dNdS.med0 = dNdS.med[,c(1,3)]
tHs0 = median(dNdS.med0[dNdS.med0>1 & dNdS.med0!= Inf], na.rm=T) # .6
tHs0 = c(tHs0, tHs0^(-1))^(-1)

##### import, restructure & filter sequence variation types #####
f.pa14 = list.files(paste0(pT[3],"../PA14_sub"), "--dbSum", full.names = T)
f.pa14.nam = read.table(text = sub("[+]","@",sub("-PA14_","@PA14_",f.pa14)), sep = "@")[,2]
for(i in 1:length(f.pa14)){
  i0 = read.csv(f.pa14[i], header = T)
  i0$PAnum = f.pa14.nam[i]
  if(i>1){dbSum = rbind(dbSum,i0)}else{dbSum = i0}
};rm(i,i0,f.pa14.nam,f.pa14)
dbSum$clinical = read.table(text = sub("_ASM","@",dbSum$clinical), sep = "@")[,1]
dbSum = dbSum[which(dbSum$clinical %in% mEta$assemblyInfo.genbankAssmAccession[mEta$sOurce %in% c("Cystic fibrosis", "Environmental")]),]

##### record & categorize median dN/dS #####
r0 = data.frame(gene = unique(dbSum$PAnum), dNdS.median.CF=NA, dNdS.median.env=NA)
for(i in 1:nrow(r0)){
  r0$dNdS.median.CF[i] = median(dbSum$dNdS[dbSum$PAnum==r0$gene[i] & dbSum$clinical %in% mEta$assemblyInfo.genbankAssmAccession[mEta$sOurce %in% "Cystic fibrosis"]], na.rm = T)
  r0$dNdS.median.env[i] = median(dbSum$dNdS[dbSum$PAnum==r0$gene[i] & dbSum$clinical %in% mEta$assemblyInfo.genbankAssmAccession[mEta$sOurce %in% "Environmental"]], na.rm = T)
};rm(i)
r0$category.CF = ifelse(is.na(r0$dNdS.median.CF),NA,ifelse(r0$dNdS.median.CF<tHs0[1],-1,ifelse(r0$dNdS.median.CF>tHs0[2],1,0)))
r0$category.env = ifelse(is.na(r0$dNdS.median.env),NA,ifelse(r0$dNdS.median.env<tHs0[1],-1,ifelse(r0$dNdS.median.env>tHs0[2],1,0)))
