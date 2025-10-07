#!/bin/env Rscript
# author: ph-u
# script: corPslMucA_stats.r
# desc: correlation between psl and other genes
# in: Rscript corPslMucA_stats.r
# out: res/corPslMucA_stats.*
# arg: 0
# date: 20251004

cat(date(),": corPslMucA_stats.r starts\n")
suppressMessages(library(Biostrings))
source("p_src.r")
d0 = read.csv(paste0(pT[1],"pslIntactness--data.csv"), header = T)[,-2] # T = psl intact

##### Record by isolates the nucleotides in enriched positions #####
iGene0 = c("PA0763", "PA0905", "PA1097", "PA1430", "PA3622")
ntLst = strsplit(nT$nucleotide[which(nT$code=="n")], ";")[[1]]

for(iG in 1:length(iGene0)){
  cat(date(),":",iGene0[iG],",",iG,"/",length(iGene0),"\n")
## identify enriched nucleotides
  sC = read.csv(paste0(pT[2],"pslIntactness--",GTF$gNam[grep(paste0("^",iGene0[iG],"$"),GTF$Locus.Tag)],"_data.csv"), header = T)
  sC = sC[!is.na(sC$oUtlier),]

  mucA = as.character(read.FASTA(paste0(pT[3],f$fNam[which(f$gEne==iGene0[iG])]), type = "DNA"))
  names(mucA)[-1] = read.table(text = gsub("_ASM", "@", names(mucA)[-1]), sep = "@")[,1]

  d1 = as.data.frame(matrix(0, nrow = nrow(d0), ncol = nrow(sC)*length(ntLst)))
  colnames(d1) = paste0(iGene0[iG],".",rep(row.names(sC), each = length(ntLst)),".",ntLst)
  row.names(d1) = d0$acc

## extract enriched nucleotides from each isolate
  for(i in 1:nrow(d1)){
    cat(date(),": checking isolates",i,"/",nrow(d1),"(",round(i/nrow(d1)*100,2),"% )","     \r")
    if(length(mucA[[which(names(mucA)==row.names(d1)[i])]]>0)){
      s0 = pairwiseAlignment(paste0(mucA[[1]], collapse = ""),paste0(mucA[[which(names(mucA)==row.names(d1)[i])]], collapse = ""), scoreOnly=F)
      s0 = as.data.frame(strsplit(as.character(c(alignedPattern(s0), alignedSubject(s0))), ""))
      if(length(grep("-",s0[,1]))>0){s0 = s0[which(s0[,1]!="-"),]}
      s0 = s0[as.numeric(row.names(sC)),]
      d1[i,paste0(iGene0[iG],".",row.names(sC),".",s0[,2])] = 1
  }};cat("\n");rm(i)
  if(iG>1){eNrich = cbind(eNrich,d1[,which(colSums(d1)>0)])}else{eNrich = d1[,which(colSums(d1)>0)]}
};rm(iG)
save(eNrich, file=paste0(pT[1],"corPslMucA_stats--data.rda"))

#### Correlation on enriched nucleotides and psl intactness #####
r0.c = c("position","estimate","rho","p.val")
r0 = as.data.frame(matrix(0, nrow = ncol(eNrich), ncol = length(r0.c)))
colnames(r0) = r0.c
r0$position = colnames(eNrich)
for(i in 1:ncol(eNrich)){
  s0 = suppressWarnings(cor.test(as.numeric(!d0$psl),eNrich[,i], method = "spearman")) # damaged psl = T
  r0[i,-1] = c(s0$estimate, s0$statistic, s0$p.value)
};rm(i, s0)
r0$p.adj = p.adjust(r0$p.val, method = "bonferroni") # bonferroni/BH
write.csv(r0[which(r0$p.adj<.05),], paste0(pT[2],"corPslMucA_stats--spearman.csv"), row.names = F, quote = F)

cat(date(),": corPslMucA_stats.r done\n")
