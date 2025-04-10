#!/bin/env Rscript
# author: ph-u
# script: src.r
# desc: env setting for whole genome dN/dS analyses
# in: source("src.r")
# out: NA
# arg: 0
# date: 20240309

library(ape)

cBp = c();for(i in c("Okabe-Ito", "alphabet", "polychrome 36", "dark 2", "set 1", "classic tableau")){cBp = c(cBp, rev(palette.colors(palette = i, alpha=1, recycle = F)))};rm(i);cBp = unique(cBp);cBp = cBp[-11]

uniProt = read.table(paste0(pT[4],"uniprotkb_taxonomy_id_208964_2024_07_05.tsv"), header = T, sep = "\t", quote = "", comment = "")
stringDB.c = strsplit(strsplit(readLines(paste0(pT[4],"208964.protein.links.full.v12.0.txt"), n=1), "[.]")[[1]], " ")[[1]]
string.DB = read.table(paste0(pT[4],"208964.protein.links.full.v12.0.txt"), header = F, quote = "", sep = " ", skip = 1)
colnames(string.DB) = stringDB.c; rm(stringDB.c)

protList = unname(unlist(read.csv(paste0(pT[1],"PAO1_proteins.txt"), header = F)))
f = list.files(pT[3],"_db.fa")
f = data.frame(fNam=f, gEne=read.table(text=gsub("[.]csv","",f), sep = "_")[,4])

library(pcaMethods) # https://doi.org/10.1093/bioinformatics/btm069, 10.18129/B9.bioc.pcaMethods
# https://stats.stackexchange.com/questions/35561/imputation-of-missing-values-for-pca
library(ggbiplot) # library(mixOmics)

##### f: capitalise first letter #####
capFirst = function(x){return(paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x))))}

##### f: use prcomp framework, ppca content in ggbiplot #####
# https://stackoverflow.com/questions/49641896/r-how-to-use-ggbiplot-with-pcares-object-plot-pca-results-of-data-with-missing
pcaSwap = function(df){
    pcaM = pca(df, method = "ppca")
    df[is.na(df)] = 999
    fRame = prcomp(df) # scale. = T : unit variance; == data normalization
    fRame$x <- pcaM@scores
    fRame$rotation <- pcaM@loadings
    fRame$sdev <- pcaM@sDev
    fRame$center <- pcaM@center
    fRame$scale <- pcaM@scale
    return(fRame)}

##### f: PCA label change #####
pcaLAB = function(pCa,nUm){
    x0 = strsplit(gsub("%","(",pCa$labels$x),"[(]")[[1]]
    pCa$labels$x = paste0(x0[1],"(",nUm[1],"%",x0[3])
    x0 = strsplit(gsub("%","(",pCa$labels$y),"[(]")[[1]]
    pCa$labels$y = paste0(x0[1],"(",nUm[2],"%",x0[3])
    return(pCa)}

##### CDS metadata #####
GTF = read.table(paste0(pT[4],"features.txt"), sep = "\t", header = T, quote = "", comment.char = "")
GTF = cbind(GTF[,-ncol(GTF)], read.table(text = sub("[.][.]","@",sub("[(]","@",GTF$Coordinates)),sep = "@")[,-3])
colnames(GTF)[(ncol(GTF)-1):ncol(GTF)] = c("start","end")
GTF$gNam0 = ifelse(GTF$Gene.Name=="",GTF$Locus.Tag,GTF$Gene.Name)

##### Mbp positions #####
m0 = floor(log10(max(GTF$end)))
Mbp = data.frame(locusTag=NA, loc=paste(1:m0,"Mbp"))
for(i in 1:m0){
    i0 = which(GTF$start>(i*10^m0))[1]
    i1 = abs(c(GTF$start[i0],GTF$end[i0-1])-i*10^m0)
    Mbp$locusTag[i] = ifelse(which(i1==min(i1))==1,GTF$NCBI.Locus.Tag[i0],GTF$NCBI.Locus.Tag[i0-1])
};rm(i,i0,i1,m0)

##### operon map #####
opMap = read.table(paste0(pT[4],"opMap--raw.gff3"), sep = "\t", quote = "")[-1,-c(6,8)]
colnames(opMap) = c("loc","src","type","start","end","strand","details")
opMap$nAm = read.table(text = sub(";","@",read.table(text = sub(";Name=","@",gsub("\'","",opMap$details)), sep = "@")[,2]), sep = "@")[,1]
opMap = opMap[grep("PA2125", opMap$nAm):grep("PA2384", opMap$nAm),]

##### STRING DB PAO1 gene name #####
string.gNam = rbind(read.table(paste0(pT[3],"../PAnumList-1.tsv"), header = T, sep = "\t", comment.char = "", quote = ""), rbind(read.table(paste0(pT[3],"../PAnumList-2.tsv"), header = T, sep = "\t", comment.char = "", quote = ""), read.table(paste0(pT[3],"../PAnumList-3.tsv"), header = T, sep = "\t", comment.char = "", quote = "")))
GTF$gNam1 = string.gNam$preferredName[match(GTF$Locus.Tag,string.gNam$queryItem)]
GTF$gNam1[is.na(GTF$gNam1)] = GTF$gNam0[is.na(GTF$gNam1)]
GTF$gNam = ifelse(GTF$gNam0==GTF$gNam1,GTF$gNam0,ifelse(substr(GTF$gNam0,1,2)=="PA",GTF$gNam1,ifelse(substr(GTF$gNam1,1,2)=="PA",GTF$gNam0,paste0(GTF$gNam0,"/",GTF$gNam1))))
GTF$gNam0 = GTF$gNam1 = NULL
rm(string.gNam)
