#!/bin/env Rscript
# author: ph-u
# script: sum_dNdS.r
# desc: statistical test for dN/dS portion of IPCD isolates
# in: Rscript sum_dNdS.r [src file name tag] [strain] [gene]
# out: data/sDNDS_[strain]_[gene].csv
# arg: 0
# date: 20240311

argv = (commandArgs(T))
# argv = c("00_PAO1_107_PA0001.csv", "PAO1_107", "PA0001")
source("metaPrep.r")
idNA = read.csv(paste0("../data/",gsub("^00_","01_",argv[1])), header = T)
id00 = idNA
id00$dNdS[grep("ID%",id00$locus)] = id00$dNdS[which(!grepl("%",idNA$locus) & is.na(idNA$dNdS))] = 0 # set 0 the isolates and gene segments that showed identical to the ref seq

##### set median collection dataframe #####
iDx0 = c("idNA","id00")
cLinc = unique(mEta$assemblyInfo.genbankAssmAccession)
dNdS.m = as.data.frame(matrix(nr = length(cLinc), nc = length(iDx0)))
row.names(dNdS.m) = cLinc ; colnames(dNdS.m) = iDx0
# dNdS.m$sRc = mEta$sOurce[match(mEta$assemblyInfo.genbankAssmAccession, row.names(dNdS.m))] # map source

##### get median of respective isolates #####
#cat("Summarizing dN/dS calculations:", date(), "\n")
for(i in 1:nrow(dNdS.m)){
#    cat(i,"/",nrow(dNdS.m),"(",round(i/nrow(dNdS.m)*100,2),"% ):",date(),"\r")
    i0 = grep(row.names(dNdS.m)[i], idNA$clinical)
    dNdS.m[i,] = c(median(idNA$dNdS[i0], na.rm = T), median(id00$dNdS[i0], na.rm = T))
};rm(i)
write.csv(dNdS.m, paste0("../data/sDNDS_",argv[2],"_",argv[3],".csv"), row.names = T, quote = F)
