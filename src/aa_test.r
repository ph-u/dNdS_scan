#!/bin/env Rscript
# author: ph-u
# script: aa_test.r
# desc: test PA2934 amino acid conservation
# in: Rscript aa_test.r
# out: NA
# arg: 0
# date: 20240730

source("p_src.r")
aa0 = as.character(read.FASTA(paste0(pT[3],f$fNam[which(f$gEne=="PA2934")]), type = "DNA"))

aa1 = as.data.frame(matrix(nr = length(aa0), nc = (length(aa0[[1]])/3)-1))
row.names(aa1) = names(aa0)
i0 = c();for(i in 1:length(aa0)){if(length(aa0[[i]])>0){aa1[i,] = strsplit(nt2prot(paste0(aa0[[i]], collapse = "")), "")[[1]]}else{i0 = c(i0,i)}};rm(i)
aa1 = aa1[-i0,];rm(i0)

aa2 = vector(mode = "list", length = nrow(aa1))
names(aa2) = row.names(aa1)
for(i in 1:nrow(aa1)){aa2[[i]] = as.character(aa1[i,])};rm(i)
write.FASTA(as.AAbin(aa2), paste0(pT[1],gsub("db","dbAA",f$fNam[which(f$gEne=="PA2934")])))
