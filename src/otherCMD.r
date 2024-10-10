#!/bin/env Rscript
# author: ph-u
# script: otherCMD.r
# desc: other commands in data exploration
# in: Rscript: otherCMD.r
# out: NA
# arg: 0
# date: 20240716

source("p_src.r")

##### Number of PAO1-specific genes #####
load(paste0(pT[1],"varTypeBar_circle.rda"))
for(i in 1:length(f.r0)){
    f.r0[[i]]$sRc = names(f.r0)[i]
    f.r0[[i]]$gEne = row.names(f.r0[[i]])
    if(i>1){r0 = rbind(r0, f.r0[[i]])}else{r0 = f.r0[[i]]}
};rm(i)

tHs = 75

cat("\nDetails on the extensive noHit fraction in the PAO1 genome:\n")
r1 = r0[which(r0$noHit>=tHs),c("gEne", "sRc","noHit")]
r2 = as.data.frame(matrix(nr = length(unique(r1$gEne)), nc = length(unique(r1$sRc))))
colnames(r2) = unique(r1$sRc)
row.names(r2) = unique(r1$gEne)
for(i in 1:nrow(r1)){ r2[which(row.names(r2) == r1$gEne[i]),which(colnames(r2) == r1$sRc[i])] = r1$noHit[i] };rm(i)
r2 = round(r2,1)
r2 = r2[which(!is.na(rowSums(r2[,1:2]))),1:2]
#r2[is.na(r2)] = paste0("< ",tHs)
cat("There were",nrow(r2),"genes\n")
print(r2)

cat("\nDetails on the extensive indel fraction in the PAO1 genome:\n")
r0$INDEL.sel = r0$INDEL.identical + r0$INDEL.SNP
r1 = r0[which(r0$INDEL.sel>=tHs),c("gEne", "sRc","INDEL.sel")]
r3 = as.data.frame(matrix(nr = length(unique(r1$gEne)), nc = length(unique(r1$sRc))))
colnames(r3) = unique(r1$sRc)
row.names(r3) = unique(r1$gEne)
for(i in 1:nrow(r1)){ r3[which(row.names(r3) == r1$gEne[i]),which(colnames(r3) == r1$sRc[i])] = r1$INDEL.sel[i] };rm(i)
r3 = round(r3,1)
r3 = r3[which(!is.na(rowSums(r3[,1:2]))),1:2]
cat("There were",nrow(r3),"genes\n")
print(r3)

##### Lost-gene island #####
print(wilcox.test(f.r0[[1]][which(f.r0[[1]]$gEne=="PA2125"):which(f.r0[[1]]$gEne=="PA2384"),"noHit"], f.r0[[2]][which(f.r0[[2]]$gEne=="PA2125"):which(f.r0[[2]]$gEne=="PA2384"),"noHit"]))

##### Neutral-selected genes #####
cat("\nDetails on neutral-selection along the PAO1 genome, CF vs Environment only:\n")
dNdS.cir = read.csv(paste0(pT[1],"dNdS_median-cir.csv"), header = T)
dC.0 = dNdS.cir[which(dNdS.cir$dNdS==2),-3]
dC.table = as.data.frame.matrix(table(dC.0))
dC.table = dC.table[which(rowSums(dC.table)<ncol(dC.table)),c(1,3)]
dC.table = dC.table[which(rowSums(dC.table)==1),]
cat("\nList of genes neutrally-selected only in CF (",nrow(dC.table[dC.table[,1]==1,]),"):",paste0(paste0(row.names(dC.table[dC.table[,1]==1,]), sep = ", "), collapse = ""),"\n")
cat("\nList of genes neutrally-selected only in the environment (",nrow(dC.table[dC.table[,2]==1,]),"):",paste0(paste0(row.names(dC.table[dC.table[,2]==1,]), sep = ", "), collapse = ""),"\n")

##### Count core genome members #####
# http://www.pnas.org/cgi/doi/10.1073/pnas.1419677112
eSs = read.table(paste0(pT[4],"pnas.1419677112.sd01.csv"), sep = "\t", header = T)
eSs = eSs[,c(1,grep("Essen", colnames(eSs)))]
eSs[eSs=="Essential"] = 0 
eSs[eSs=="Nonessential"] = 1 
for(i in 2:ncol(eSs)){eSs[,i] = as.numeric(eSs[,i])};rm(i)

eSs.0 = rowSums(eSs[,-1])
eSs.c = c("condition","gene.count",paste0(rep(c("c","e"), each = 3),".",c("cons", "neut", "vari"), sep = ""))
eSs.1 = as.data.frame(matrix(nr = 4, nc = length(eSs.c)))
colnames(eSs.1) = eSs.c
eSs.1[,1] = c("Hard core", "Conditionally-core", "Soft core", "Non-essential")
eSs.1[,2] = c(sum(eSs.0==0), sum(eSs.0>0 & eSs.0<(ncol(eSs)-1)), sum(eSs.0<(ncol(eSs)-1)), sum(eSs.0==(ncol(eSs)-1)))

dC.src = c("1.Cystic.fibrosis", "3.Environmental")
gType = list(eSs$Locus.ID[which(eSs.0==0)], eSs$Locus.ID[which(eSs.0>0 & eSs.0<(ncol(eSs)-1))], eSs$Locus.ID[which(eSs.0<(ncol(eSs)-1))], eSs$Locus.ID[which(eSs.0==(ncol(eSs)-1))]) # list of hard core, Conditionally-core, soft core, non-essential
dC.ce = dNdS.cir[which(dNdS.cir$sRc %in% dC.src),]
for(i in 1:nrow(eSs.1)){ for(i1 in 1:length(dC.src)){
    i4 = dC.ce[which(dC.ce$sRc==dC.src[i1] & dC.ce$gEne %in% gType[[i]]),]
    i0 = table(i4$dNdS)
    i2 = rep(0,length(unique(dC.ce$dNdS[!is.na(dC.ce$dNdS)])))
    i2[as.numeric(names(i0))] = i0
    eSs.1[i,2+length(i2)*(i1-1)+(1:length(i2))] = i2
    if(i<3){
        i4 = i4[which(i4$dNdS > 1),]; i4 = i4[order(i4$dNdS),-2]
        i4$dNdS = ifelse(i4$dNdS>2, "vari", "neut")
        cat(eSs.1[i,1],"in",dC.src[i1],"\n")
        print(i4)
}}};rm(i, i0, i1, i2, i4)
eSs.1$c.abs = eSs.1[,2]-(eSs.1[,3]+eSs.1[,4]+eSs.1[,5])
eSs.1$e.abs = eSs.1[,2]-(eSs.1[,6]+eSs.1[,7]+eSs.1[,8])
eSs.1 = eSs.1[,c(1,2,9,3:5,10,6:8)]
print(eSs.1)

cat("Essential genes: CF vs env\n")
print(chisq.test(as.numeric(eSs.1[3,3:6]/eSs.1[3,2]), as.numeric(eSs.1[3,7:10]/eSs.1[3,2])))
cat("Non-essential genes: CF vs env\n")
print(chisq.test(as.numeric(eSs.1[4,3:6]/eSs.1[4,2]), as.numeric(eSs.1[4,7:10]/eSs.1[4,2])))
cat("CF: Essential vs non\n")
print(chisq.test(as.numeric(eSs.1[3,3:6]/eSs.1[3,2]), as.numeric(eSs.1[4,3:6]/eSs.1[4,2])))
cat("env: Essential vs non\n")
print(chisq.test(as.numeric(eSs.1[3,7:10]/eSs.1[3,2]), as.numeric(eSs.1[4,7:10]/eSs.1[4,2])))

# dNdS.cir[which(is.na(dNdS.cir$dNdS) & dNdS.cir$sRc==dC.src[1] & dNdS.cir$gEne %in% gType[[1]]),] # show individual dNdS type content
# intersect(dNdS.cir[which(dNdS.cir$dNdS==1 & dNdS.cir$sRc==dC.src[2] & dNdS.cir$gEne %in% gType[[1]]),1], dNdS.cir[which(dNdS.cir$dNdS==3 & dNdS.cir$sRc==dC.src[1] & dNdS.cir$gEne %in% gType[[1]]),1]) # show comparative difference between two gene categories

# https://doi.org/10.1073/pnas.1900570116
#eS0.nam = paste0(pT[1],"pnas.1900570116.modData.csv")
#if(!file.exists(eS0.nam)){
#    eS0 = read.csv(paste0(pT[4],"pnas.1900570116.sd05.csv"), header = T)
#    eS.eq = read.csv(paste0(pT[4],"Pseudomonas_aeruginosa_UCBPP-PA14_109_orthologs.csv"), header = T)
#    eS.eq = eS.eq[grep("PAO1", eS.eq$Strain.Hit.),]
#    eS0$PAnum = eS.eq$Locus.Tag.Hit.[match(eS0$PA14_ID, eS.eq$Locus.Tag.Query.)]
#    write.csv(eS0, eS0.nam, row.names = F, quote = F)
#    rm(eS.eq, eS0)
#};eS0 = read.csv(eS0.nam, header = T)
#eS0 = eS0[!is.na(eS0$PAnum.QC),]
#eSs$Essential..3 = ifelse(eS0$Essential.Category[match(eSs$Locus.ID, eS0$PAnum.QC)]=="Core",0,1)
#eSs$Essential..3[is.na(eSs$Essential..3)] = 1

##### Correlation between high mutational burden & dN/dS #####
vOl = read.csv(paste0(pT[1],"dNdS_infect-vol.csv"), header = T)
vOl = vOl[which(vOl$cond1=="Cystic fibrosis" & vOl$cond2=="Environmental"),]
v.0 = rep(NA, 3)

## being essential correlate with significant CF-env dN/dS discrepancy?
cat("Hard core genes (CF/env ratio) vs whole dN/dS background:\n")
print(v.1 <- wilcox.test(vOl$log2FC[which(vOl$gene %in% gType[[1]])], vOl$log2FC)) # hard core
v.0[1] = v.1$p.value

cat("Soft core genes (CF/env ratio) vs whole dN/dS background:\n")
print(v.1 <- wilcox.test(vOl$log2FC[which(vOl$gene %in% gType[[3]])], vOl$log2FC)) # soft core
v.0[2] = v.1$p.value

## being high mutational burden in CF correlate with significant CF-env dN/dS discrepancy?
cat("High mutational burden genes (CF/env ratio) vs whole dN/dS background:\n")
hL = paste0("PA",c("0424","0425","0762","0763","0958","0861","1288","1430","1798",2020,2426,2491,2492,3168,3458,3477,3545,3974,4266,4367,4522,4601,5291), sep = "")
print(v.1 <- wilcox.test(vOl$log2FC[which(vOl$gene %in% hL)], vOl$log2FC))
v.0[3] = v.1$p.value
print(p.adjust(v.0, method = "bonf"))
rm(v.0, v.1)

##### List of genes conserved in CF and environment respectively #####
cat("\nList of genes more conserved in CF strains:\n")
paste0(vOl$gene[which(vOl$log2FC < -1 & vOl$p.adj < .1)], collapse = ",")
cat("\nList of genes more conserved in environmental strains:\n")
paste0(vOl$gene[which(vOl$log2FC > 1 & vOl$p.adj < .1)], collapse = ",")
