#!/bin/env Rscript
# author: ph-u
# script: cif_phy.r
# desc: phylogenetic analysis on PA2934 (Cif) contrasting Pseudomonas ortholog group POG000798
# in: Rscript cif_phy.r
# out: NA
# arg: 0
# date: 20240813

source("p_src.r")
library(Biostrings);library(msa)
orthGp = paste0(pT[4],"POG000798",".",c("fasta","tab"), sep = "")

##### Combine fasta files #####
cif.r = as.character(read.FASTA(orthGp[1], type = "AA"))
cif.o = as.character(read.FASTA(paste0(pT[4],"alphaBetaHydrolase.fasta"), type = "AA"))
for(i in 1:length(cif.o)){cif.r[length(cif.r)+1] = cif.o[i];names(cif.r)[length(cif.r)] = names(cif.o)[i]};rm(i)
write.FASTA(as.AAbin(cif.r), paste0(pT[1],"cif_phy.fa"))

# https://doi.org/10.1371/journal.pone.0243927
# https://www.r-bloggers.com/2012/07/r-for-ecologists-phylogenies-in-r/
# http://stackoverflow.com/questions/21565143/ddg#21566066
cat("Multiple sequence alignment:",date(),"\n")
cif.r = msa(readAAStringSet(paste0(pT[1],"cif_phy.fa")), method = "ClustalOmega", verbose = T)
cat("Conversion to AAbin:",date(),"\n")
cif.aabin = as.AAbin(cif.r)
cat("Constructing phylogeny:",date(),"\n")
cif.phy = nj(dist.aa(cif.aabin))
cif.phy = root(cif.phy, grep("Luteimonas abyssi",cif.phy$tip.label))

cat("Writing phylogeny:",date(),"\n")
cif.phy0 = cif.phy
cif.phy0$tip.label = paste0(cif.phy0$tip.label,"_",1:length(cif.phy0$tip.label))
write.tree(cif.phy0, paste0(pT[1],"cif_phy.nwk"))

cat("Drawing phylogeny:",date(),"\n")
#cif.phy = cif.phy0
otgData = read.table(orthGp[2], header = T, sep = "\t")
otgData$ePi = NA; for(i in 1:nrow(otgData)){otgData$ePi[i] = strsplit(otgData$Strain[i], " ")[[1]][2]};rm(i)
pPat = read.csv(paste0(pT[4],"pseudomonasPathogenicity.csv"), header = T)
otgData$pAt = pPat$pathogenReported[match(otgData$ePi,pPat$Epithet)]; otgData$pAt[is.na(otgData$pAt)] = "Unknown"
otgData$ePi0 = 1+as.numeric(as.factor(otgData$ePi))
otgData$pAt0 = 1+as.numeric(as.factor(otgData$pAt))
#for(i in grep("Assembly", otgData$Strain)){otgData$Strain[i] = strsplit(sub(" - Assembly", "@", otgData$Strain[i]), "@")[[1]][1]};rm(i)
#for(i in grep("contigs", otgData$Strain)){otgData$Strain[i] = strsplit(sub("_contigs", "@", otgData$Strain[i]), "@")[[1]][1]};rm(i)
#otgData$nAm = paste(gsub("Pseudomonas ","P. ",otgData$Strain),otgData$Locus.Tag,"ln",row.names(otgData))
#cif.tips = read.table(text = sub("[|]","@",read.table(text = gsub("[|]gi[|]","@",cif.phy$tip.label), sep = "@")[,2]), sep = "@")[,1]
#otgData$sEq = NA; for(i in 1:nrow(otgData)){otgData$sEq[i] = grep(otgData$Locus.Tag[i],cif.phy$tip.label)[1]};rm(i)
#cif.phy$tip.label[-grep("Altererythrobacter|Streptomyces",cif.phy$tip.label)] = otgData$nAm[otgData$sEq]

#pdf(paste0(pT[2],"cif_phy-tre.pdf"), height = 500, width = 50)
#plot(cif.phy)
#invisible(dev.off())

cif.o.nam = paste0(read.table(text = sub("[]]","",names(cif.o)), sep = "[")[,2],collapse = "|")
pdf(paste0(pT[2],"cif_phy-cir.pdf"), height = 40, width = 30)
par(mfrow = c(2,1), mar = c(0,0,0,0)+.1)
for(i0 in 1:2){if(i0==1){
        cOl = rep(1,length(cif.phy$tip.label)); for(i in 1:nrow(otgData)){cOl[grep(otgData$RefSeq.Accession[i], cif.phy$tip.label)] = otgData$ePi0[i]};rm(i)
    }else{
        cOl = rep(1,length(cif.phy$tip.label)); for(i in 1:nrow(otgData)){cOl[grep(otgData$RefSeq.Accession[i], cif.phy$tip.label)] = otgData$pAt0[i]};rm(i)
    }
# https://rdrr.io/cran/ape/man/plot.phylo.html
    plot.phylo(cif.phy, type = "fan", label.offset = 1, no.margin = T, cex = .1, tip.color = cBp[cOl])
    for(i in grep("PA0829|PA2934|PA2949|PA3053|PA3429|PA3994|PA4152",cif.phy$tip.label)){tiplabels(strsplit(cif.phy$tip.label[i]," ")[[1]][1], i, frame = "n", col = "#ff00ffff", cex = 2)};rm(i)
    for(i in grep(cif.o.nam,cif.phy$tip.label)){tiplabels(sub("[]]","",strsplit(cif.phy$tip.label[i], "[[]")[[1]][2]), i, frame = "n", col = ifelse(i==grep("vinaceusdrappus",cif.phy$tip.label), "#0055ffff", "#0000ffff"), cex = .5)};rm(i)
};rm(i0)
invisible(dev.off())
#write.tree(cif.phy, paste0(pT[1],"cif_phy-tre.nwk"))

pdf(paste0(pT[2],"cif_phy-cirL.pdf"), height = 14, width = 7)
par(mfrow = c(1,2), mar = c(0,0,0,0)+.1)
cif.L = unique(otgData[,c("ePi","ePi0")])
plot(x = c(0,0), y = c(0,mM<-nrow(cif.L)), col = "#00000000")
text(x = 0, y = 1:nrow(cif.L), labels = cif.L$ePi, col = cBp[cif.L$ePi0])
cif.L = unique(otgData[,c("pAt","pAt0")])
plot(x = c(0,0), y = c(0,max(mM,nrow(cif.L))), col = "#00000000")
text(x = 0, y = 1:nrow(cif.L), labels = cif.L$pAt, col = cBp[cif.L$pAt0])
invisible(dev.off())

cat("cif_phy.r completed:",date(),"\n")

# https://askubuntu.com/questions/1011336/does-there-exist-a-pdf-viewer-for-ubuntu-offering-up-to-6400-percent-zoom#1083653
# gsettings set org.gnome.Evince page-cache-size 5000
