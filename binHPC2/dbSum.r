#!/bin/env Rscript
# author: ph-u
# script: dbSum.r
# desc: summarize database "_db.fa" blastn searches on one gene
# in: Rscript dbSum.r [inFile_db.fa]
# out: data/[strain]_[gene]--dbSum.csv
# arg: 1
# date: 20240330 (branch from sumBLASTN.sh)

argv = (commandArgs(T))
#argv = "1_GCA_000022165_STM14_0013"
if(length(grep("_db.fa", argv))>0){inFile=argv;setwd("/binHPC2")}else{inFile=0} # docker
cat(argv,":",date(),"\n")

source("src_dNdS.r"); suppressMessages(library(ape)); suppressMessages(library(Biostrings))
pT = paste0("../",c("data","binHPC2"),"/")
aRg = list.files(pT[1],argv)
argv = aRg[grep("_db.fa",aRg)]

##### set work env #####
## fixed variables
db.flk = 100 # flanking region length

## import
cat(date(),": dbSum.r data import start\n")
if(inFile==0){
  f = as.character(read.FASTA(paste0(pT[1], argv[1]), type="DNA"))
  fNam = sub("_db.fa","",argv)
  fNam = strsplit(fNam,"--")[[1]]
  fNam = c(fNam, sub(paste0("_",fNam[2]),"",fNam[1]))
  oNam = paste0(pT[1], paste0(fNam[3:2], collapse = "--"), "--dbSum.csv")
  lGenome=".."
}else{ # docker
  argv = inFile
  f = as.character(read.FASTA(argv, type="DNA"))
  fNam = strsplit(sub("-","+",sub("_db.fa","",sub("../data/","",argv))),"[+]")[[1]]
  oNam = sub("_db.fa","--dbSum.csv",argv)
  lGenome="../data"
}
listGenome = list.files(lGenome, recursive=T, include.dirs=F, full.names=T)
if(inFile==0){
  fLank = read.csv(paste0(pT[1],fNam[3],"-flanking.csv"), header = F)
  pao1 = which(fLank$PAnum==fNam[2])
}else{
  fLank = read.csv(listGenome[grep("flank",listGenome)][grep(fNam[1],listGenome[grep("flank",listGenome)])], header = F)
  pao1 = which(fLank$PAnum==fNam[2] & paste0(fLank$start,fLank$end, sep = "")==fNam[3])
}
colnames(fLank) = c("index","PAnum","start","end","beforeGFlank","afterGFlank")

## PAO1 details
pao1 = list(seq=paste0(f[[1]], collapse = ""), details=c(length(f[[1]]), as.numeric(fLank[pao1,3:4])), flanks=as.character(fLank[pao1,5:6]))

## set rec df
cat(date(),": dbSum.r Setting data record format\n")
r0.c = c("clinical","locus","varType","start","end","flank_befPc","flank_aftPc","dNdS","pN","pS","Nd","Sd","N","S")
r0 = as.data.frame(matrix(nr = length(f)-1, nc = length(r0.c)))
colnames(r0) = r0.c
if(file.exists(oNam)){
    r0.p = read.csv(oNam, header = T)
    i0 = nrow(r0.p)+1
    r0[1:nrow(r0.p),] = r0.p
    rm(r0.p)
}else{i0=1}

##### Process each db #####
cat(date(),": dbSum.r Processing each database\n")
if(i0 <= nrow(r0)){ for(i in i0:nrow(r0)){
    db.nam = strsplit(names(f)[i+1], ";")[[1]]
    db.nRec = strsplit(sub("_geno","@",db.nam[1]), "@")[[1]][1]
    cat(date(),": ",i,"/",nrow(r0),"(",round(i/nrow(r0)*100),"% ;",db.nRec,")", date(), "\r")
    if(db.nam[2]=="noHit"){ r0[i,c(1,3)] = c(db.nRec,db.nam[2]) }else{
        r0[i,1:2] = c(db.nRec,db.nam[2])

## Sequence variation type identification
        r0[i,4:5] = db.stEd = as.numeric(strsplit(db.nam[3], "-")[[1]])
        db.seq = paste0(f[[i+1]], collapse = "")
	if(pao1$seq==db.seq){r0$varType[i] = "identical"}else{
            db.sep = pairwiseAlignment(pao1$seq, db.seq)
            db.sep = c(as.character(db.sep@pattern), as.character(db.sep@subject))
            if(length(grep("-",db.sep))==0 & nchar(db.sep[1])==nchar(db.sep[2])){
                r0$varType[i] = "SNP"
                r0[i,grep("dNdS", colnames(r0)):ncol(r0)] = dNdS.rt(db.sep[1], db.sep[2])[[1]]
            }else{
## Identify indel type: count long seq insert -> codon -> 2bp -> 1bp
                vt.frameshift = sum(strsplit(paste0(db.sep, collapse = ""), "")[[1]]=="-") %% 3
                vt.fs4 = sum(unlist(strsplit(vt.seq <- gsub("----+","@",db.sep),"")) == "@")
                vt.fs3 = sum(unlist(strsplit(vt.seq <- gsub("---","@",gsub("@","",vt.seq)),"")) == "@")
                vt.fs2 = sum(unlist(strsplit(vt.seq <- gsub("--","@",gsub("@","",vt.seq)),"")) == "@")
                vt.fs1 = sum(unlist(strsplit(gsub("-","@",gsub("@","",vt.seq)),"")) == "@")
                vt.t1 = "indel"
## Consider if indel does not interfere with ORF or more than 90% sequence intect (C-terminus frameshifted)
# https://stackoverflow.com/questions/11619616/how-to-split-a-string-into-substrings-of-a-given-length
                if(vt.frameshift == 0 | (length(grep("-",db.sep))==0 & nchar(db.sep[1])*.9 <= nchar(db.sep[2]))){
                    vt.df = strsplit(db.sep, "") # 2-steps seq->codons
                    vt.df = data.frame(pao1=paste0(vt.df[[1]][c(T,F,F)],vt.df[[1]][c(F,T,F)],vt.df[[1]][c(F,F,T)]), db=paste0(vt.df[[2]][c(T,F,F)],vt.df[[2]][c(F,T,F)],vt.df[[2]][c(F,F,T)]))
                    vt.df = vt.df[-unique(c(grep("-",as.data.frame(t(vt.df))),which(nchar(vt.df[,1]) %% 3 != 0),which(nchar(vt.df[,2]) %% 3 != 0))),]
### calcaulate dN/dS with approx indel location(s) removed
                    db.sep = c(paste0(vt.df[,1], collapse = ""), paste0(vt.df[,2], collapse = ""))
		    if(db.sep[1]==db.sep[2]){ vt.t1 = "identical" }else{ vt.t1 = "SNP"
                        r0[i,grep("dNdS", colnames(r0)):ncol(r0)] = dNdS.rt(db.sep[1], db.sep[2])[[1]]
		    }};r0$varType[i] = paste0("INDEL:FS",vt.frameshift,":",vt.fs1,":",vt.fs2,":",vt.fs3,":",vt.fs4,":",vt.t1) } }

## flanking region similarity
        db.stEd = db.stEd + c(ifelse(db.nam[4]=="plus",-1,1),ifelse(db.nam[4]=="plus",1,-1))
        db.stEd = c(db.stEd[1]+(db.flk-1)*ifelse(db.nam[4]=="plus",-1,1), db.stEd[1], db.stEd[2], db.stEd[2]+(db.flk-1)*ifelse(db.nam[4]=="plus",1,-1))

        db.f = as.character(read.FASTA(listGenome[grep(db.nam[1],listGenome)], type="DNA"))
        db.f = db.f[[grep(db.nam[2], names(db.f))]]

## db flanking regions integrity check
        db.stEd[db.stEd<1] = 1
        db.stEd[db.stEd>length(db.f)] = length(db.f)

## extract flanking regions & calculate similarity score (default scoring ref from Biostrings)
# https://evomics.org/resources/substitution-models/nucleotide-substitution-models/
        if(db.stEd[2]>1 & db.stEd[2]<length(db.f)){
            db.seq = paste0(db.f[db.stEd[1]:db.stEd[2]], collapse = "")
            r0$flank_befPc[i] = pairwiseAlignment(tolower(substr(pao1$flanks[1], db.flk-nchar(db.seq)+1, nchar(pao1$flanks[1]))),tolower(db.seq), gapOpening = 0, gapExtension = 0, scoreOnly=T)/pairwiseAlignment(db.seq,db.seq, gapOpening = 0, gapExtension = 0, scoreOnly=T)*100 }
        if(db.stEd[3]>1 & db.stEd[3]<length(db.f)){
            db.seq = paste0(db.f[db.stEd[3]:db.stEd[4]], collapse = "")
            r0$flank_aftPc[i] = pairwiseAlignment(tolower(substr(pao1$flanks[2], db.flk-nchar(db.seq)+1, nchar(pao1$flanks[2]))),tolower(db.seq), gapOpening = 0, gapExtension = 0, scoreOnly=T)/pairwiseAlignment(db.seq,db.seq, gapOpening = 0, gapExtension = 0, scoreOnly=T)*100 }

##### export #####
    };write.csv(r0[!is.na(r0[,1]),], oNam, row.names = F, quote = F) } }
cat(argv[1], "overall summary done:", date(), "\n")
