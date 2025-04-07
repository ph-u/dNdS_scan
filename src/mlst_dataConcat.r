#!/bin/env Rscript
# author: ph-u
# script: mlst_dataConcat.r
# desc: Concatenate MLST sequences for phylogeny building
# in: Rscript mlst_dataConcat.r {1..7}
# out: [data_folder]/mlst/*.fna
# arg: 0
# date: 20250407

argv = (commandArgs(T))
source("p_src.r")
oNam = paste0(pT[3],"../ncbi-genomes-2023-05-15/mlst/mlst--",argv[1],".fa")
library(msa) # https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
f0 = f[which(f$gEne %in% paste0("PA",c("0025", "0609", "0887", "1770", "2639", "3769", "4946"))),]

i=as.numeric(argv[1])
cat(date(),": Start alignment, ",i,"\n")
#for(i in 1:nrow(f0)){
  cat(date(),": Aligning ",f0$gEne[i],", ",i,"/",nrow(f0),"\n")
  i0 = as.character(read.FASTA(paste0(pT[3],f0$fNam[i]), type = "DNA"))
#  i0 = i0[1:5] # testing
  if(length(grep("noHit", names(i0)))>0){for(i1 in grep("noHit", names(i0))){i0[[i1]] = rep("n",length(i0[[1]]))}}
  names(i0) = c("PAO1", read.table(text = sub("_ASM","@",names(i0)[-1]), sep = "@")[,1])
  write.FASTA(as.DNAbin(i0), oNam)
  r0 = as.list(as.character(msa(readDNAStringSet(oNam), "ClustalW")))
#  if(i>1){
#    for(i1 in 1:length(r0)){r0[[i1]] = paste0(r0[i1],aLn[which(names(aLn)==names(r0)[i1])])}
#  }else{r0 = aLn}
#};rm(i,i0,i1)
cat("\n",date(),": Start exporting, ",i,"\n")
r00 = strsplit(as.character(r0),"")
names(r00) = names(r0)
#write.FASTA(as.DNAbin(r00), paste0(pT[1],"mlst_raw.fa"))
write.FASTA(as.DNAbin(r00), oNam)
#if (file.exists(oNam)) {file.remove(oNam)}
cat(date(),": MLST raw data exported, ",i,"\n")
