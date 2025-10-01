#!/bin/env Rscript
# author: ph-u
# script: p_corPslMucA.r
# desc: Does damage of psl operon intactness correlate with mucA conservation?
# in: Rscript p_corPslMucA.r
# out: res/p_corPslMucA.tif
# arg: 0
# date: 20250723

set.seed(20250724)
library(msa)
source("p_src.r")
paLst = paste0("PA",c(2231:2245,"0763"))

cat(date(),": sequence variation data extraction\n")
for(i in 1:length(paLst)){
  cat(date(),": fusing data     \r")
  i0 = read.csv(paste0(pT[3],"PAO1_107_",paLst[i],"--dbSum.csv"), header = T)[,c(1,3)]
  if(length(grep("INDEL",i0[,2]))>0){
    i0[grep("INDEL",i0[,2]),2] = apply(read.table(text = i0[grep("INDEL",i0[,2]),2], sep = ":")[,c(1,7)], 1, function(s){return(paste0(s[1],":",s[2]))})
  }
  colnames(i0)[2] = GTF$gNam[which(GTF$Locus.Tag==paLst[i])]
#  i0[,2] = as.factor(i0[,2])
  if(i>1){r0 = merge(r0,i0)}else{r0 = i0} 
};rm(i,i0);cat("\n",date(),": fusing completed\n")

##### Extracting clinical isolates #####
r0$clinical = read.table(text = gsub("_ASM", "@", r0$clinical), sep = "@")[,1]
r0.c = r0[which(r0$clinical %in% mEta$assemblyInfo.genbankAssmAccession[grep("Cystic",mEta$sOurce)]),]
r0.cIntact = data.frame(psl = apply(r0.c[,2:(ncol(r0.c)-1)], 1, function(x){length(grep("noHit", x))==0}), mucA = (r0.c$mucA %in% c("SNP","INDEL:SNP")))

##### mucA & lasR gene between intact and non-intact psl operon #####
iGene0 = c("PA0763", "PA1430")
for(iG in length(iGene0)){
  mucA = as.character(read.FASTA(paste0(pT[3],f$fNam[which(f$gEne==iGene0[iG])]), type = "DNA"))
  names(mucA)[-1] = read.table(text = gsub("_ASM", "@", names(mucA)[-1]), sep = "@")[,1]

  seq.pslIntact = mucA[which(names(mucA) %in% r0.c$clinical[which(r0.cIntact$psl)])]
  seq.pslLoss = mucA[which(names(mucA) %in% r0.c$clinical[which(!r0.cIntact$psl)])]

  cat("Isolates that loss the",iGene0[iG],"\nIntact psl operon: ")
  i0 = c();for(i in 1:length(seq.pslIntact)){if(length(seq.pslIntact[[i]])==0){i0 = c(i0,i);cat(names(seq.pslIntact)[i],",")}};rm(i);cat(";; Total number =",length(i0),"\nDamaged psl operon: ")
  if(length(i0)>0){seq.pslIntact = seq.pslIntact[-i0]}
  i0 = c();for(i in 1:length(seq.pslLoss)){if(length(seq.pslLoss[[i]])==0){i0 = c(i0,i);cat(names(seq.pslLoss)[i],",")}};rm(i);cat(";; Total number =",length(i0),"\n")
  if(length(i0)>0){seq.pslLoss = seq.pslLoss[-i0]}

  write.FASTA(as.DNAbin(seq.pslIntact), paste0(pT[1],iGene0[iG],"-pslIntact.fa"))
  write.FASTA(as.DNAbin(seq.pslLoss), paste0(pT[1],iGene0[iG],"-pslLoss.fa"))

## Seq alignment between two groups
  cat(date(),": Multiple Sequence Alignment - psl partial loss in CF-associated strains\n")
  seq.pslLoss = msa(readDNAStringSet(paste0(pT[1],iGene0[iG],"-pslLoss.fa")), method = "ClustalW")
  cat(date(),": Multiple Sequence Alignment - psl intact in CF-associated strains\n")
  seq.pslIntact = msa(readDNAStringSet(paste0(pT[1],iGene0[iG],"-pslIntact.fa")), method = "ClustalW")
  cat(date(),"MSA done\n")

  cQ.L = msaConsensusSequence(seq.pslLoss)
  cQ.I = msaConsensusSequence(seq.pslIntact)

  write.FASTA(as.DNAbin(msaConvert(seq.pslIntact, type = "seqinr::alignment")), paste0(pT[1],iGene0[iG],"-pslIntact.fa"))
  write.FASTA(as.DNAbin(msaConvert(seq.pslLoss, type = "seqinr::alignment")), paste0(pT[1],iGene0[iG],"-pslLoss.fa"))

## Plot nucleotide conservation
  sL = as.data.frame(as.character(read.FASTA(paste0(pT[1],iGene0[iG],"-pslLoss.fa"), type = "DNA")))
  sI = as.data.frame(as.character(read.FASTA(paste0(pT[1],iGene0[iG],"-pslIntact.fa"), type = "DNA")))

  seq.c = c("a","t","c","g")
  sL.df = as.data.frame(matrix(0, nrow = nrow(sL), ncol = length(seq.c)))
  sI.df = as.data.frame(matrix(0, nrow = nrow(sI), ncol = length(seq.c)))
  colnames(sL.df) = colnames(sI.df) = seq.c

  cat(date(),": Counting started\n")
  for(i in 1:max(nrow(sL),nrow(sI))){cat(date(),": ",i,"/",max(nrow(sL),nrow(sI)),"   \r")
    if(i<=nrow(sL)){for(i0 in 1:length(seq.c)){sL.df[i,i0] = sum(sL[i,]==seq.c[i0])}}
    if(i<=nrow(sI)){for(i0 in 1:length(seq.c)){sI.df[i,i0] = sum(sI[i,]==seq.c[i0])}}
  };rm(i);cat("\n",date(),": done\n")

  sL.df = sL.df/ncol(sL)
  sI.df = sI.df/ncol(sI)

  sL.df$seq = strsplit(cQ.L, "")[[1]]
  sI.df$seq = strsplit(cQ.I, "")[[1]]

  sI.df = sI.df[which(sI.df$seq!="-"),]
  sL.df = sL.df[which(sL.df$seq!="-"),]

  seqContrast = sI.df[, -ncol(sI.df)]/sL.df[, -ncol(sL.df)]

# default 1:1 ratio
  seqContrast[seqContrast==0] = 1
  seqContrast[is.na(seqContrast)] = 1
  for(i in 1:ncol(seqContrast)){seqContrast[is.infinite(seqContrast[,i]),i] = 1};rm(i)
  seqContrast = log10(seqContrast)

  tHs = .1
  seqContrast$oUtlier = NA
  for(i in 1:nrow(seqContrast)){if(any(seqContrast[i,-ncol(seqContrast)] > tHs)){seqContrast$oUtlier[i] = .3}else if(any(seqContrast[i,-ncol(seqContrast)] < -tHs)){seqContrast$oUtlier[i] = -.3}}

  tHs = which(!is.na(seqContrast$oUtlier))
  jpeg(paste0(pT[2],ifelse(iG==1,"mucA","lasR"),"_pslIntactness.jpeg"), width = 2000, height = 1500, res = 300)
  par(mar = c(5,4,0,1)+.1)
  matplot(x = 1:nrow(seqContrast), y = seqContrast[,-ncol(seqContrast)], type = "l", lty = 1, lwd = 2, col = cBp, xlab = "Sequence order", ylab = "log10( Ratio of occurrence )")
  text(x = tHs, y = seqContrast$oUtlier[tHs]*(runif(length(tHs))+.5), labels = tHs)
  text(x = c(50,50), y = .3*c(-1,1), labels = c("With\ngene\nloss", "Intact"), cex = 2)
  legend(ifelse(iG==1,"topright","bottomleft"), legend = capFirst(seq.c), fill = cBp)
  invisible(dev.off())

##### Mutations in mucA #####
## N/S type?
  mucA.mut = mucA[[1]]
  site.mut = seqContrast[tHs,]
  for(i in 1:nrow(site.mut)){
    mucA.mut[tHs[i]] = seq.c[which.max(abs(site.mut[i,-ncol(site.mut)]))]
  };rm(i)
  print(mucA.dNdS <- dNdS.rt(paste0(mucA[[1]], collapse = ""), paste0(mucA.mut, collapse = "")))

  write.csv(mucA.dNdS[[2]], paste0(pT[2],ifelse(iG==1,"mucA","lasR"),"_pslIntactness.csv"), row.names = F, quote = F)
};rm(iG, iGene0)
