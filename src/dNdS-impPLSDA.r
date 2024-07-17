#!/bin/env Rscript
# suthor: ph-u
# script: dNdS-impPLSDA.r
# desc: Imputated dN/dS data using PLSDA
# in: Rscript dNdS-impPLSDA.r [num_of_NA_allowed]
# out: {data,res}/dNdS-impPLSDA.*
# arg: 1
# date: 20240624

argv = (commandArgs(T))
tRs = ifelse(argv[1]=="", 9, as.numeric(argv[1]))
source("p_src.r")
library(missForest)
f.rda = "dNdS-impPLSDA.rda"
if(tRs>nrow(mEta)){stop(paste0("Max threshold = ",nrow(mEta),"\n"))}

if(!file.exists(paste0(pT[1],f.rda))){
f0 = list.files(pT[3], "dbSum")
f = data.frame(fNam = f0, gEne = read.table(text = gsub("_PA","@PA",gsub("--","@",f0)), sep = "@")[,2])
cat("Start extracting dN/dS data into one data.frame:", date(), "\n")
for(i in 1:nrow(f)){ cat(i,"/",nrow(f),"(",round(i/nrow(f)*100,2),"% )", date(), "          \r")
    i0 = read.csv(paste0(pT[3],f$fNam[i]), header = T)
    if(i==1){
        r0 = as.data.frame(matrix(nr=nrow(i0), nc=nrow(f)+1))
        row.names(r0) = i0$clinical
        colnames(r0) = c("sRc",f$gEne)
        i0$clinical = read.table(text = gsub("_ASM","@",i0$clinical), sep = "@")[,1]
        r0[,1] = mEta$sOurce[match(i0$clinical, mEta$assemblyInfo.genbankAssmAccession)]
        r0[,i+1] = i0$dNdS
    }else{
        r0[,i+1] = i0$dNdS[match(row.names(r0), i0$clinical)]
}};rm(i, i0);cat("\nDone collection:", date(),"\n")
save(r0, file = paste0(pT[1],f.rda))
}else{load(paste0(pT[1],f.rda))}

##### Exclude genes with all NAs #####
if(!file.exists(paste0(pT[1],gsub("-imp",paste0("_",tRs,"-imp"),f.rda)))){
for(i in 2:ncol(r0)){r0[is.infinite(r0[,i]),i] = NA};rm(i) # treat dN/dS not applicable data as missing
x0 = rep(NA,nrow(f))
cat("Get genes with NAs:",date(),"\n")
for(i in 1:length(x0)){cat(i,":",date(),"     \r");x0[i] = sum(is.na(r0[,i+1]))};rm(i);cat("\nDone scanning NA genes:", date(), "\n")
r0 = r0[,c(1,which(x0<tRs)+1)] # exclude amount of missing data above threshold
cat("Take",sum(x0<tRs), "genes (out of ",length(x0),"), NA allowance threshold = ",tRs,":",date(),"\n")

##### Plot NA histogram #####
cat("Plot NA count gene frequency:",date(),"\n")
x00 = table(x0)
pdf(paste0(pT[2], gsub("rda","pdf",gsub("-impPLSDA",paste0("_",tRs,"-impHist"),f.rda))), width = 14, height = 7)
plot(x=names(x00), y = x00, main="Distribution of Missing data and\nNon-synonymous only Genes", xlab = "NA count in genes", ylab = "Gene Frequency", pch = 3)
abline(v=tRs, col = cBp[1], lty = 2, lwd = 3)
invisible(dev.off())

##### Imputation missing dN/dS values #####
cat("Imputation start:",date(),"\n")
r0.0 = missForest(r0[,-1])$ximp
cat("Imputation done:",date(),"\n")
r0.0$sRc = r0$sRc
r0.0 = r0.0[,c(ncol(r0.0),1:(ncol(r0.0)-1))]
save(r0.0, file = paste0(pT[1],gsub("-imp",paste0("_",tRs,"-imp"),f.rda)))
}else{load(paste0(pT[1],gsub("-imp",paste0("_",tRs,"-imp"),f.rda)))}

##### PLS-DA #####
library(mixOmics)
r.plsda = plsda(r0.0[,-1], as.factor(r0.0$sRc))
pdf(paste0(pT[2],gsub("rda","pdf",gsub("-imp",paste0("_",tRs,"-imp"),f.rda))), width = 6, height = 7)
plotIndiv(r.plsda, pch=20, ellipse=T, ellipse.level=.95, col = cBp[1:length(unique(r0.0$sRc))], size.xlabel=rel(1.4), size.ylabel=rel(4), size.axis=rel(2), legend=T, style="graphics", cex=1)
invisible(dev.off())

# for i in `seq 11 14`;do Rscript dNdS-impPLSDA.r ${i} >> ../data/00_dNdS_${i}-impPLSDA.txt;done
