#!/bin/env Rscript
# author: ph-u
# script: clinical_clusLost.r
# desc: map Lost genes with STRING cluster
# in: Rscript clinical_clusLost.r
# out: NA
# arg: 0
# date: 20240823, 20250211

source("p_src.r")
library(lattice);library(ggplot2);library(stringr)

load(paste0(pT[1],"varTypeBar_circle.rda"))
for(i in 1:length(f.r0)){
    f.r0[[i]]$sRc = names(f.r0)[i]
    f.r0[[i]]$gEne = row.names(f.r0[[i]])
    if(i>1){r0 = rbind(r0, f.r0[[i]])}else{r0 = f.r0[[i]]}
};rm(i)

lGene = read.csv(paste0(pT[1],"clinical_gMap.csv"), header = T, row.names = 1)
colnames(lGene) = read.table(text = gsub("_ASM","@",colnames(lGene)), sep = "@")[,1]
i.0 = mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==unique(mEta$sOurce)[grep("Cystic",unique(mEta$sOurce))])]
stringClus = read.table(paste0(pT[4],"string_kmeans_clusters--PA.tsv"), sep = "\t", header = T, comment.char = "", quote = "")

##### data trim: CF & within lost gene region #####
lGene = lGene[which(row.names(lGene) %in% stringClus$protein.name),which(colnames(lGene) %in% i.0)]
lG0 = lGene<1 # gene lost
print(stringClus[which(stringClus$protein.name %in% row.names(lG0)[which(rowSums(lG0)>(.7*ncol(lG0)))]),c(2,3,6)])

for(i in 1:nrow(lG0)){
    lG0[i,] = ifelse(lG0[i,],stringClus$cluster.number[which(stringClus$protein.name==row.names(lG0)[i])],0)
};rm(i)

pdf(paste0(pT[2],"clinical_clusLost--heat.pdf"), width = 10)
print(levelplot(t(lG0), col.regions = rep(c("#00000000",paste0(unique(stringClus$hex.color),"77", sep = "")), each = 4), xlab = "clinical isolates from individuals with Cystic Fibrosis", ylab = "PAO1 Gene Number (PA2125 -> PA2384)")) # cBp[c(3,5,6,7,4)]
invisible(dev.off())

##### operon map #####
opMap$yy = 0
repeat{ p0 = nrow(opMap)
    for(i in 2:nrow(opMap)){opMap$yy[i] = ifelse(opMap$end[i-1]>opMap$start[i],opMap$yy[i-1]+1,opMap$yy[i])};rm(i)
#    opMap = opMap[opMap$yy == 0,]
    if(nrow(opMap)==p0){rm(p0);break}
}
i0 = strsplit(opMap$nAm, "-")
opMap$nAmPlot = ""
for(i in 1:nrow(opMap)){if(length(grep("PA",i0[[i]]))!=length(i0[[i]]) | length(grep("PA",i0[[i]]))==0){opMap$nAmPlot[i] = opMap$nAm[i]}};rm(i,i0)

opMap$yyTxt = 0
opMap$yyTxt[opMap$nAmPlot!=""] = rep(c(1,-1),sum(opMap$nAmPlot!=""))[1:sum(opMap$nAmPlot!="")]

##### Patches within contingency island #####
lG1 = data.frame(PAnum = row.names(lG0), gNam = GTF$gNam[match(row.names(lG0), GTF$Locus.Tag)], lostRatio = rowSums(lG0>0)/ncol(lG0), cluster = stringClus$cluster.color[match(row.names(lG0), stringClus$protein.name)])
lG1$operon = NA
for(i in 1:nrow(lG1)){if(length(grep(lG1$gNam[i],opMap$nAm))>0){lG1$operon[i] = opMap$nAm[grep(lG1$gNam[i],opMap$nAm)[1]]}};rm(i)

pdf(paste0(pT[2],"clinical_clusLost--gLostRatio.pdf"), width = 10, height = 7)
par(mar = c(5,5,0,0)+.1, cex.axis = 3, cex.lab = 3)
plot(x = 1:nrow(lG1),y = lG1$lostRatio*100, xaxt = "n", pch = 20, col = lG1$cluster, xlab = "PA2125 -> PA2384", ylab = paste0("gene lost (%, n = ",ncol(lG0),")"))
abline(h = c(summary(f.r0[[1]][-(which(f.r0[[1]]$gEne == "PA2125"):which(f.r0[[1]]$gEne == "PA2384")),"noHit"])["Mean"], summary(f.r0[[2]][which(f.r0[[2]]$gEne=="PA2125"):which(f.r0[[2]]$gEne=="PA2384"),"noHit"])["Mean"]))
invisible(dev.off())
print("High percentage lost CDS")
print(lG1[which(lG1$lostRatio>.3),])

##### plot operon map #####
opMap$color = lG1$cluster[match(opMap$nAm, lG1$operon)]
opMap$color[is.na(opMap$color)] = c("Red","Red","Yellow","Yellow","Yellow","Cyan","Cyan")
lGCol = unique(stringClus[,c("cluster.color","hex.color")])
lGCol[nrow(lGCol)+1,] = c(NA,"#000000ff")
p0 = ggplot() + theme_minimal() + xlab("") + ylab("") + ylim(c(-1,5)) + coord_radial(r.axis.inside = T, inner.radius = .7) +
     geom_segment(aes(x = opMap$start[opMap$strand=="+"], y=opMap$yy[opMap$strand=="+"], xend = opMap$end[opMap$strand=="+"], yend = opMap$yy[opMap$strand=="+"], colour = opMap$color[opMap$strand=="+"]), arrow = arrow(length = unit(.2, "cm"), type = "closed"), show.legend = F) +
     geom_segment(aes(x = opMap$end[opMap$strand!="+"], y=opMap$yy[opMap$strand!="+"], xend = opMap$start[opMap$strand!="+"], yend = opMap$yy[opMap$strand!="+"], colour = opMap$color[opMap$strand!="+"]), arrow = arrow(length = unit(.2, "cm"), type = "closed"), show.legend = F) +
     scale_colour_manual(values = setNames(object = lGCol$hex.color,lGCol$cluster.color)) +
#     geom_segment(aes(x = rowSums(opMap[opMap$nAmPlot!="",c("start","end")])/2, y=opMap$yy[opMap$nAmPlot!=""], xend = rowSums(opMap[opMap$nAmPlot!="",c("start","end")])/2, yend = opMap$yyTxt[opMap$nAmPlot!=""]*.5), arrow = arrow(length = unit(.2, "cm"))) +
#     geom_text(aes(x = rowSums(opMap[,c("start","end")])/2, y = opMap$yyTxt, label = opMap$nAmPlot), angle = 0) +
     theme(axis.text = element_blank())

ggsave(paste0(pT[2],"clinical_clusLost--operonMap.pdf"), plot = p0, width = 10, height = 10)

