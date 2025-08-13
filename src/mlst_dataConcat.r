#!/bin/env Rscript
# author: ph-u
# script: mlst_dataConcat.r
# desc: Concatenate MLST sequences for phylogeny building
# in: Rscript mlst_dataConcat.r
# out: [data_folder]/mlst/mlst--all.fa
# arg: 0
# date: 20250408

source("p_src.r")
oNam = paste0(pT[3],"../ncbi-genomes-2023-05-15/mlst/MLST--all.fa")

rUn = 0
##### Data Concat #####
if(rUn>0){
f.mlst = list.files(paste0(pT[3],"../ncbi-genomes-2023-05-15/mlst"),"mlst--", full.names = T)
for(i in 1:length(f.mlst)){
  i0 = as.character(read.FASTA(f.mlst[i], type = "DNA"))
  if(i>1){
    i1 = match(names(r0),names(i0)) # mapping i0 to r0
    for(i2 in 1:length(r0)){r0[[i2]] = c(r0[[i2]],i0[[i1[i2]]])}
  }else{r0 = i0}
};rm(i,i0,i1,i2)
write.FASTA(as.DNAbin(r0), oNam)
write.FASTA(as.DNAbin(r0[-which(names(r0)=="PAO1")]), sub("all","IPCD",oNam))
write.FASTA(as.DNAbin(r0[which(names(r0) %in% mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce %in% c("Cystic fibrosis", "Environmental"))])]), sub("all","cfEnv",oNam))
}

##### Phylogeny #####
#pHy.all = njs(dist.dna(read.FASTA(sub("all","IPCD",oNam)), model = "TN93"))
#pHy.CF = njs(dist.dna(read.FASTA(sub("all","cfEnv",oNam)), model = "TN93"))
a = as.character(read.FASTA(sub("all","cfEnv",oNam), type = "DNA"))
write.FASTA(as.DNAbin(a[-grep("GCA_015277575.1",names(a))]), sub("all","cfEnvnoOut",oNam));rm(a)
pHy.CF = njs(dist.dna(read.FASTA(sub("all","cfEnvnoOut",oNam)), model = "TN93"))
#tipCol = rep("#000000FF",length(pHy.CF$tip.label))
#tipCol[which(pHy.CF$edge.length==max(pHy.CF$edge.length))+(-4:-1)] = cBp[3]
#tipCol[1:4] = cBp[3]

#highEvoClade = paste0("GCA_003",c(83987,83902,83724,84065,83482,83975,83719,97388,83764,83702,97508,83670,84020,83368,97381,83652),"5.1")
#hEc.meta = mEta[which(mEta$assemblyInfo.genbankAssmAccession %in% pHy.CF$tip.label[c(1:4,which(pHy.CF$edge.length==max(pHy.CF$edge.length))+(-4:-1))]),c(3,16,17,14,15)]
#for(i in ncol(hEc.meta)-(0:1)){hEc.meta[,i] = gsub(",",";",hEc.meta[,i])};rm(i)
#write.csv(hEc.meta, paste0(pT[2],"mlst--highEvo.csv"), quote = F, row.names = F)

#pdf(paste0(pT[2],"mlst--phylogeny.pdf"), width = 12, height = 24, paper = "a4")
#par(mfrow = c(4,1))
#plot(pHy.CF, "f", main = paste0("MLST tree of CF- / Environmental-associated isolates; n = ",length(pHy.CF$tip.label)), cex=.1, tip.color = tipCol)
#tipCol = cBp[as.numeric(as.factor(mEta$sOurce[match(pHy.CF$tip.label,mEta$assemblyInfo.genbankAssmAccession)]))]
#plot(pHy.CF, "f", main = paste0("MLST tree of CF- / Environmental-associated isolates; n = ",length(pHy.CF$tip.label)), cex=.1, tip.color = tipCol)
#legend("bottomleft", legend = levels(as.factor(mEta$sOurce[match(pHy.CF$tip.label,mEta$assemblyInfo.genbankAssmAccession)])), fill = cBp[1:length(unique(tipCol))])
#tipCol = cBp[as.numeric(as.factor(mEta$cOuntry[match(pHy.CF$tip.label,mEta$assemblyInfo.genbankAssmAccession)]))]
#plot(pHy.CF, "f", main = paste0("MLST tree of CF- / Environmental-associated isolates; n = ",length(pHy.CF$tip.label)), cex=.1, tip.color = tipCol)
#legend("bottomleft", legend = levels(as.factor(mEta$cOuntry[match(pHy.CF$tip.label,mEta$assemblyInfo.genbankAssmAccession)])), fill = cBp[1:length(unique(tipCol))], ncol = 5)
#plot(pHy.all, "f", main = paste0("MLST tree of all isolates from the initial IPCD collection; n = ",length(pHy.all$tip.label)), cex=.5)
#invisible(dev.off())

#jpeg(paste0(pT[2],"mlst--phyCFENV.jpeg"), width = 1000, height = 300, res = 500)
#par(mai=c(0,0,0,0)+.1)
#pHy.CF0 = pHy.CF
#pHy.CF0$edge.length[which(pHy.CF0$edge.length==max(pHy.CF0$edge.length))] = max(pHy.CF0$edge.length)/5
#plot(pHy.CF0, "f", cex=.01)
#invisible(dev.off())

#jpeg(paste0(pT[2],"mlst--phyIPCD.jpeg"), width = 1000, height = 700, res = 300)
#par(mai=c(0,0,0,0)+.1)
#plot(pHy.all, "f", cex=.1)
#invisible(dev.off())

for(i in 1:2){
  jpeg(paste0(pT[2],"mlst--phyCFENV_",ifelse(i==1,"src","geo"),".jpeg"), width = 2000, height = 2000, res = 500)
  par(mai=c(0,0,0,0)+.1)
  tipCol = cBp[ifelse(i==1,2,0)+as.numeric(as.factor(mEta[match(pHy.CF$tip.label,mEta$assemblyInfo.genbankAssmAccession),ifelse(i==1,"sOurce","cOuntry")]))]
  plot(pHy.CF, "f", cex=.1, tip.color = tipCol)
  invisible(dev.off())

  jpeg(paste0(pT[2],"mlst--phyCFENV_",ifelse(i==1,"src","geo"),"L.jpeg"), width = 2000, height = 2000, res = 500)
  par(mai=c(0,0,0,0)+.1)
  plot(0,0, col = "#00000000")
  legend("topleft", legend = levels(as.factor(mEta[match(pHy.CF$tip.label,mEta$assemblyInfo.genbankAssmAccession),ifelse(i==1,"sOurce","cOuntry")])), fill = cBp[(1:length(unique(tipCol)))+ifelse(i==1,2,0)], ncol = ifelse(i==1,1,2))
  invisible(dev.off())
};rm(i,tipCol)
