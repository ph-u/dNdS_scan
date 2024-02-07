#!/bin/env Rscript
# author: ph-u
# script: rollingdNdS.r
# desc: rolling sequence window on all available Pseudomonas aeruginosa PA01 genes against one clinical sample
# in: Rscript rollingdNdS.r [#clinical]
# out: res/01_dNdS_[clinical].r
# arg: 1
# date: 20231217

argv = (commandArgs(T))
#argv=c("PAO1_107","PGD108156")
library(ape); source("src_dNdS.r")
pT = "../data/" #paste0("../",c("data","dNdSGenome"),"/")

PAfa = list.files(pT[1],"fa")
cLfa = list.files(pT[1],"_db.fa$")

pA = as.character(read.FASTA(paste0(pT[1],PAfa[grep(paste0(argv[1],"_",argv[2],".fa"), PAfa)])))[[1]]
bSum = read.csv(paste0(pT[1],"00_",argv[1],"_",argv[2],".csv"), header = T)

##### Result collection preparation #####
aaLen = 67
aaLen = ifelse(length(pA) > (aaLen+1)*3,aaLen,ifelse(length(pA) > (aaLen-aaLen%%2)/2*3,(aaLen-aaLen%%2)/2,ceiling(aaLen/10))) # 3-layered filter: 67, 33, 7 amino acids length window
i = nrow(bSum)
x = c("gene", "clinical", "locus", "start", "end", "flank_befPc", "flank_aftPc", "dNdS","pN","pS","Nd","Sd","N","S")

if(length(grep("%",bSum$locus[i]))==0 & length(pA) > (aaLen+1)*3){
    cSeq = as.character(read.FASTA(paste0(pT[1],cLfa[grep(paste0(argv[1],"_",argv[2]), cLfa)])))
    c0 = cSeq[[grep(bSum$clinical[i], names(cSeq))[1]]]
    if(bSum$dbSeqDir[i]=="minus"){for(i0 in 1:length(c0)){c0[i0] = nFlip(c0[i0])}}

##### dNdS calculation #####
    cat("Calculating dNdS for:",bSum[i,1],"against",bSum[i,2],date(),"\n")

## construct data-recording dataframe
    stPos = seq(1, min(length(pA),length(c0))-(aaLen-1)*3, 3)
    dRec = as.data.frame(matrix(nr = length(stPos), nc = length(x)))
    dRec[,1:3] = bSum[i,1:3]
    dRec[,6:7] = bSum[i,7:8]
    dRec[,4] = stPos
    dRec[,5] = dRec[,4] + aaLen*3-1
    dRec[which(dRec[,5] > length(pA) | dRec[,5] > length(c0)),5] = min(length(pA),length(c0))

## segment sequences
    for(i1 in 1:nrow(dRec)){
        a1 = paste0(pA[dRec[i1,4]:dRec[i1,5]], collapse = "") # PA01
	if(bSum$dbSeqDir[i]=="minus"){# clinical
	    a2 = paste0(rev(c0[dRec[i1,4]:dRec[i1,5]]), collapse = "") # minus
	}else{
            a2 = paste0(c0[dRec[i1,4]:dRec[i1,5]], collapse = "") # plus
	}
        dRec[i1,-c(1:(grep("dNdS",x)-1))] = dNdS.rt(a1, a2)[[1]]
    }; rm(i1, a1, a2)
}else{
    dRec = as.data.frame(matrix(nr = 1, nc = length(x)))
    dRec[,1:4] = bSum[i,c(1:3,5)]
    dRec[,3] = ifelse(length(pA) > (aaLen+1)*3,bSum[i,3],paste0("SHORT%",bSum[i,3]))
    dRec[,5] = dRec[,4]+bSum[i,6]
    dRec[,6:7] = bSum[i,7:8]
}
cat("Exporting:",date(),"\n")
write.csv(dRec, paste0(pT[1], "01_",argv[1],"_",argv[2],"_t.csv"), quote = F, row.names = F)
cat("Completed:",date(),"\n")
