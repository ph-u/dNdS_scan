#!/bin/env Rscript
# author: ph-u
# script: rDNDS.r
# desc: rolling sequence window dN/dS calculations on one gene
# in: Rscript rDNDS.r [inFile_db.fa]
# out: data/[strain]_[gene]--rDNDS.csv
# arg: 1
# date: 20240331 (branch from rollingdNdS.r)

argv = (commandArgs(T))
#argv = "00_PAO1_107_PA0001_db.fa"
if(length(grep("_db.fa", argv))>0){inFile=argv;setwd("/binHPC2")}else{inFile=0} # docker
cat(argv,":",date(),"\n")

source("src_dNdS.r"); suppressMessages(library(ape)); suppressMessages(library(Biostrings))
pT = paste0("../",c("data","binHPC2"),"/")
aRg = list.files(pT[1],argv)
argv = aRg[grep("_db.fa",aRg)]

##### set work env #####
## fixed variables
aaLen.df = 67

## import
cat(date(),": rDNDS.r Importing data\n")
if(inFile==0){
  f = as.character(read.FASTA(paste0(pT[1], argv[1]), type="DNA"))
  fNam = sub("_db.fa","",argv)
  fNam = strsplit(fNam,"--")[[1]]
  fNam = c(fNam, sub(paste0("_",fNam[2]),"",fNam[1]))
  oNam = paste0(pT[1], paste0(fNam[3:2], collapse = "--"), "--rDNDS.csv")
}else{ # docker
  argv = inFile
  f = as.character(read.FASTA(argv, type="DNA"))
  fNam = strsplit(sub("-","+",sub("_db.fa","",sub("../data/","",argv))),"[+]")[[1]]
  oNam = sub("_db.fa","--rDNDS.csv",argv)
}

dbSum = read.csv(gsub("rDNDS","dbSum",oNam), header = T)
rfSeq = paste0(f[[1]], collapse = "")

## set rec start point
cat(date(),": rDNDS.r Setting data record format\n")
r0.c = c("clinical","locus","segVarType","ntStart","ntEnd","dNdS")
if(file.exists(oNam) & file.info(oNam)$size>0){
    r0.p = read.csv(oNam, header = T)
    if(!setequal(colnames(r0.p),r0.c)){i0=1; wRite = 0}else{
        i0 = which(paste0(dbSum$clinical,dbSum$locus)==paste0(r0.p$clinical[nrow(r0.p)],r0.p$locus[nrow(r0.p)]))+1
        wRite = 1
};rm(r0.p)}else{i0=1; wRite = 0}

##### Process each db #####
cat(date(),": rDNDS.r Processing sequence sliding windows\n")
if(i0 <= nrow(dbSum)){ for(i in i0:nrow(dbSum)){ cat(i,"/",nrow(dbSum),"(",round(i/nrow(dbSum)*100),"% ;", dbSum$clinical[i],")", date(),"\n"); if(length(grep("SNP", dbSum$varType[i])) > 0 & dbSum$varType[i]!="noHit"){
    aaLen = aaLen.df
    db.sep = c(rfSeq, paste0(f[[i+1]], collapse = ""))

    if(dbSum$varType[i] != "SNP" | nchar(db.sep[1]) != nchar(db.sep[2])){ # realign seq if indel exist / partial blastn match
        db.se1 = pairwiseAlignment(db.sep[1], db.sep[2]) # avoid trimming at N-terminus
        db.se2 = pairwiseAlignment(db.sep[2], db.sep[1]) # avoid trimming at N-terminus
        db.sep = c(as.character(db.se1), as.character(db.se2))
    }

    vt.df = strsplit(db.sep, "") # 2-steps seq->codons
    vt.df = data.frame(ref=paste0(vt.df[[1]][c(T,F,F)],vt.df[[1]][c(F,T,F)],vt.df[[1]][c(F,F,T)]), db=paste0(vt.df[[2]][c(T,F,F)],vt.df[[2]][c(F,T,F)],vt.df[[2]][c(F,F,T)]))
    vt.df$ntPos = seq(1,nchar(db.sep[1]),3)

    if(dbSum$varType[i] != "SNP" | length(grep("-",vt.df))>0){ # remove indel regions & frameshifted remnants
        if(length(grep("-",vt.df[,1]))>0){ for(i1 in 2:nrow(vt.df)){ vt.df$ntPos[i1] = vt.df$ntPos[i1-1] + ifelse(length(grep("-",vt.df[i1,1]))>0,0,3) } } # collapse altered ref seq segments
        vt.df = vt.df[-unique(c(grep("-",as.data.frame(t(vt.df))),which(nchar(vt.df[,1]) %% 3 != 0),which(nchar(vt.df[,2]) %% 3 != 0))),]
    }

    aaLen = ifelse(nrow(vt.df) > aaLen*1.5, aaLen, ifelse(nrow(vt.df) > aaLen/2, floor(aaLen/2), ceiling(aaLen/10))) # 3-layer filter for genes of different lengths
    r0 = as.data.frame(matrix(nr = nrow(vt.df)-aaLen+1, nc = length(r0.c)))
    colnames(r0) = r0.c
    r0$clinical = dbSum$clinical[i]
    r0$locus = dbSum$locus[i]
    r0$ntStart = vt.df$ntPos[1:(nrow(vt.df)-(aaLen-1))]
    r0$ntEnd = vt.df$ntPos[aaLen:nrow(vt.df)]+2
    r0$segVarType = "identical"; r0$dNdS = 0 # set default value; an identical segment in an SNP gene is assumed dN/dS == 0
    for(i1 in 1:nrow(r0)){
        db.seq = c(paste0(vt.df[i1:(i1+aaLen-1),1], collapse = ""), paste0(vt.df[i1:(i1+aaLen-1),2], collapse = ""))
        if(db.seq[1]!=db.seq[2]){r0$segVarType[i1] = "SNP"; r0$dNdS[i1] = dNdS.rt(db.seq[1],db.seq[2])[[1]][1]}
    }

##### export #####
    cat("Exporting",i,"/",nrow(dbSum),"(",dbSum$clinical[i],")",date(),"\n")
    if(wRite==0){write.csv(r0, oNam, row.names = F, quote = F); wRite = 1}else{write.csv(rbind(read.csv(oNam, header = T),r0), oNam, row.names = F, quote = F)} } } }
cat(argv[1], "sequence window dN/dS calculation done:", date(), "\n")
