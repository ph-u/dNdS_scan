#!/bin/env Rscript
# author: ph-u
# script: blastn2faByGene.r
# desc: Reformat blastn result into fasta format; match blastn result with env strain genomes
# in: Rscript blastn2faByGene.r [../relative/path/blastnRES].txt [envGenome_dir]_iDx.txt
# out: [../relative/path/blastnRES]--[geneName].fa
# arg: 2|3
# date: 20241119, 20241127, 20250301

argv = (commandArgs(T))
#argv = c("../data/1_GCA_000022165_STM14_2638.txt", "1_GCA_000022165_iDx.txt")
library(ape)

blastnRES = read.table(argv[1], sep = "@", header = F)
blastnRES[,2] = read.table(text=blastnRES[,2], sep = "|")[,2]

if(length(argv)==3){ ## docker container
    srcG = as.character(read.FASTA(argv[2], type = "DNA"))
    envPT = paste0("../data/", read.table("../data/iDx.txt")[as.numeric(argv[3]),1], "/", collapse = "")
    envG = list.files(envPT, "fna")
    envG0 = list.files(envPT, "fna", full.names = T)
    rm(envPT)
}else{
    srcG = as.character(read.FASTA(paste0(dirname(argv[1]),"/",sub("_iDx.txt$",".fa",argv[2])), type = "DNA"))
    envG = list.files(read.table(argv[2])[1,1],"fna")
    envG0 = list.files(read.table(argv[2])[1,1],"fna", full.names = T)
}

##### Map contig names with genome name #####
cat("Contig-Genome mapping starts:",date(),"\n")
for(i in 1:length(envG)){
  cat(date(),":",i,"/",length(envG),"(",round(i/length(envG)*100,2),"% )     \r")
  e0 = names(as.character(read.FASTA(envG0[i], type = "DNA")))
  if(i>1){
    envNode = rbind(envNode,data.frame(fNam=sub(".fna$","",envG[i]),nNam=e0))
  }else{
    envNode = data.frame(fNam=sub(".fna$","",envG[i]),nNam=e0)
  }
};rm(i,e0);cat("\n")
cat("Contig-Genome mapping completed:",date(),"\n")

##### Map Genome name with matched sequence from blastn #####
uniQ = unique(blastnRES[,1])
dOne = c(); for(i in 1:length(uniQ)){
  br0 = blastnRES[which(blastnRES[,1] == uniQ[i]),] # get blastn result for one gene
  sg0 = grep(uniQ[i], names(srcG), fixed = T)
  dOne = c(dOne,sg0) # record which gene is completed
  gNam = strsplit(sub("]","@",sub(";","@",strsplit(sub("locus_tag=","@",names(srcG)[sg0]),"@")[[1]][2])),"@")[[1]][1]
  dna0 = strsplit(br0[,6],"") # non-ref sequence collection bin

cat("Genome matching starts:",date(),"\n")
## Match genome names
  br0$genome = apply(br0,1,function(x){return(grep(x[2], envNode$nNam))})
  names(dna0) = paste0(envNode$fNam[br0$genome],";",br0[,2],";",br0[,3],"-",br0[,4],";",br0[,5])

## Add in noHit genomes if any
  x0 = unique(envNode$fNam[br0$genome])
  if(length(x0)==length(envG)){if(all(sub(".fna$","",envG)[order(envG)]==x0[order(x0)])){e0=1}else{e0=0}}else{e0=0} # e0==0 if there are noHits
  if(e0==0){
    x0 = paste0(unique(envNode$fNam[!(envNode$fNam %in% x0)]),";noHit")
    x0 = setNames(as.list(rep("N",length(x0))),x0)
  }
cat("Genome matching done:",date(),"\n")
  if(e0==0){ r0 = c(srcG[sg0], dna0, x0) }else{ r0 = c(srcG[sg0], dna0) }
    if(length(argv)==3){ # docker container
      write.FASTA(as.DNAbin(r0), sub("[.]fa","_db.fa",sub("_","-",argv[2])))
    }else{
      write.FASTA(as.DNAbin(r0),sub(".txt$",paste0("--",gNam,"_db.fa"),argv[1]))
    }
};rm(i)
cat("blastn2faByGene.r done:",date(),"\n")
