#!/bin/env Rscript
# author: ph-u
# script: getnoHitBySrc.r
# desc: summarize gene details of the noHit dominated fraction
# in: Rscript getnoHitBySrc.r
# out: data/getnoHitBySrc.csv
# arg: 0
# date: 20240505

i0=as.numeric((commandArgs(T)))
load("../data/varTypeBar_circle.rda")
fEat = read.table("../raw/features.txt", header = T, quote = "", comment.char="", sep = "\t")

a.0 = vector(mode = "list", length = length(f.r0))
a.1 = c(); for(i in 1:length(f.r0)){
    print(nrow(a0<-f.r0[[i]][which(f.r0[[i]]$noHit>=(i0*100)),]))
    a.0[[i]] = row.names(a0)
    a.1 = unique(c(a.1,row.names(a0)))
};rm(i, a0)
a.1 = a.1[order(a.1)]
a0 = as.data.frame(matrix("", nr = length(a.1), nc = length(f.r0)))
colnames(a0) = names(f.r0)
a.1 = cbind(PAnum=a.1[order(a.1)], a0);rm(a0)

for(i in 1:length(a.0)){a.1[,i+1] = f.r0[[i]]$noHit[match(a.1$PAnum, row.names(f.r0[[i]]))]};rm(i,a.0)
#for(i in 1:length(a.0)){a.1[match(a.0[[i]],a.1$PAnum),i+1] = f.r0[[i]]$noHit[match(a.0[[i]], row.names(f.r0[[i]]))]};rm(i,a.0)
a.1$mean = apply(a.1[,-1], 1, mean)
a.1$tgt = ""
for(i in 2:nrow(a.1)){if(as.numeric(gsub("[.]1","",substr(a.1$PAnum[i-1],3,nchar(a.1$PAnum[i-1]))))+1==as.numeric(gsub("[.]1","",substr(a.1$PAnum[i],3,nchar(a.1$PAnum[i]))))){a.1$tgt[i] = "v"}};rm(i)
a.1$prod = gsub(",",";",fEat$Product.Description[match(a.1$PAnum,fEat$Locus.Tag)])

write.csv(a.1, "../data/getnoHitBySrc.csv", row.names = F, quote = F)
