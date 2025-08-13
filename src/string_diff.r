#!/bin/env Rscript
# author: ph-u
# script: string_diff.r
# desc: post hoc analysis of string k-means clustering result
# in: Rscript string_diff.r
# out: NA
# arg: 0
# date: 20240728

source("p_src.r")
library(ggdendro)
i0 = list.files(paste0(pT[1],"string_diff"), "tsv")
for(i in 1:length(i0)){
    r.r = read.table(paste0(pT[1],"string_diff/",i0[i]), sep = "\t", header = T, quote = "", comment.char = "")[,c(5,2)]
    colnames(r.r)[2] = paste0(colnames(r.r)[2],".",i)
    if(i>1){r0 = merge(r0, r.r, by = "protein.name")}else{r0 = r.r}
};rm(i, i0)
row.names(r0) = r0[,1]
r0[,1] = NULL
fEat = read.table(paste0(pT[4],"features.txt"), sep = "\t", header = T, quote = "", comment.char = "")[,c("Locus.Tag","Gene.Name")]
vOl = read.csv(paste0(pT[1],"dNdS_infect-vol.csv"), header = T)
vOl = vOl[which(vOl$cond1=="Cystic fibrosis" & vOl$cond2=="Environmental"),]

##### Hierarchical Cluster Analysis #####
hCa = hclust(dist(r0, method = "euclidean"), method = "ward.D2") # https://stats.stackexchange.com/questions/109949/what-algorithm-does-ward-d-in-hclust-implement-if-it-is-not-wards-criterion#109958

p0 = ggdendrogram(hCa, rotate = F, size = 1) +
    theme_minimal() + scale_y_reverse() + coord_radial(r.axis.inside = T, inner.radius = .3) +
    geom_text(aes(x=1:nrow(r0), y=-10, label = hCa$labels[hCa$order]), angle = ((1:nrow(r0))-1)/nrow(r0)*-360+90, size = 1) +
    theme(axis.text.x = element_blank())
ggsave(paste0(pT[2], "string_diff.pdf"), plot = p0, width = 10, height = 10)
#pdf(paste0(pT[2], "string_diff.pdf"), width = 10, height = 10)
#plot(hCa, cex = .5, hang = -1)
#invisible(dev.off())
hCut = data.frame(cID = cutree(hCa, k=20))
hCut$gEne = row.names(hCut)
hCut = hCut[order(hCut$cID),]
hCut$PAnum = ifelse(substr(hCut$gEne,1,2)=="PA",hCut$gEne,NA)
for(i in which(is.na(hCut$PAnum))){hCut$PAnum[i] = ifelse(length(which(fEat$Gene.Name==hCut$gEne[i]))==0,NA,fEat$Locus.Tag[which(fEat$Gene.Name==hCut$gEne[i])])};rm(i)

#xx0 = which(is.na(hCut$PAnum))
hCut$PAnum[is.na(hCut$PAnum)] = c("PA4790", "PA0515", "PA5254", "PA5129", "PA1197", "PA0779", "PA4006", "PA5127", "PA2606", "PA2605", "PA0503", "PA5348", "PA2963", "PA0209", "PA3388", "PA0688") # STRING DB
#print(hCut[xx0,])

hCut$CFenv = vOl$log2FC[match(hCut$PAnum, vOl$gene)]
hCut$CFenv = ifelse(hCut$CFenv < 0,"CF","env")

hCut = hCut[order(hCut[,1],hCut[,4],hCut[,3]),]
#hCut$stringDB = gsub(",","&",string.gNam$annotation[match(hCut$PAnum,string.gNam$queryItem)])
#hCut$GTF = gsub(",","&",GTF$Product.Description[match(hCut$PAnum,GTF$Locus.Tag)])
#hCut$Uniprot = gsub(",","&",uniProt$Protein.names[match(GTF$UniProtKB.ID[match(hCut$PAnum,GTF$Locus.Tag)], uniProt$Entry.Name)])
write.csv(hCut, paste0(pT[2], "string_diff.csv"), row.names = F, quote = F)
for(i in 1:max(hCut$cID)){cat("Cluster",i,"\n");print(paste0(hCut$PAnum[which(hCut$cID==i)],collapse = ","))};rm(i)
