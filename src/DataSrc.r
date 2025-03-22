#!/bin/env Rscript
# author: ph-u
# script: DataSrc.r
# desc: plot metadata country and sample sources as heatmap
# in: Rscript DataSrc.r
# out: res/DataSrc.pdf
# arg: 0
# date: 20230529

source("p_src.r")
library(lattice);library(ggmap);library(scatterpie);library(ggrepel)

x0 <- table(mEta[which(mEta$sOurce %in% c("Cystic fibrosis","Environmental")),c("cOuntry","sOurce")])
x1 = data.frame(cOuntry = row.names(x0), cOntinent = "Europe")
x1[which(x1[,1]=="Unknown"),2] = "Unknown"
x1[which(x1[,1] %in% c("Australia", "New Zealand")),2] = "Oceania"
x1[which(x1[,1] %in% c("Japan", "South Korea", "India", "Thailand", "Philippines")),2] = "Asia"
x1[which(x1[,1] %in% c("USA", "Canada", "Mexico")),2] = "America.N"
x1[which(x1[,1] %in% c("Brazil", "Colombia", "Panama", "Puerto Rico")),2] = "America.S"
x1[which(x1[,1] %in% c("Benin", "Republic of the Congo", "Tunisia")),2] = "Africa"
x1[which(x1[,1] %in% c("Georgia", "Israel")),2] = "Middle East"

##### Background distribution #####
pdf(paste0(pT[2],"DataSrc--heat.pdf"), width = 10, height = 4)
print(levelplot(x0[order(x1$cOntinent),], scales=list(x=list(rot=90)), col.regions = gray(100:0/100), ylab = "Sample source", xlab = "Country"))
invisible(dev.off())

##### Source mapping #####
map.world = map_data("world")
x0 = as.data.frame.matrix(x0)
rownames(x0)[grep("United Kingdom", rownames(x0))] = "UK"
rownames(x0)[grep("Republic of the Congo", rownames(x0))] = "Republic of Congo"
x0$radius = log10(rowSums(x0))*10
x0$radius = ifelse(x0$radius < 2,2,x0$radius)
x0$cOuntry = rownames(x0)
for(i in 1:nrow(x0)){
    x0$long[i] = mean(map.world$long[which(map.world$region==x0$cOuntry[i])])
    x0$lat[i] = mean(map.world$lat[which(map.world$region==x0$cOuntry[i])])
};rm(i)
x0[which(x0$cOuntry=="New Zealand"),"long"] = 175
x0[which(x0$cOuntry=="Canada"),"lat"] = 55
x0[which(x0$cOuntry=="Unknown"),c("long","lat")] = c(-10,-40)

cOl = cBp[1:length(unique(mEta$sOurce))]
names(cOl) = unique(mEta$sOurce)

p0 = ggplot()+theme_bw()+coord_equal()+
    geom_point(data = map.world, aes(x=long, y=lat), cex=.1, color="#00000033")+
    geom_label_repel(aes(x=long, y=lat, label=cOuntry), data=x0, size=7, max.overlaps = getOption("ggrepel.max.overlaps", default = 40))+
    geom_scatterpie(aes(x=long, y=lat, group=cOuntry, r=radius), data=x0, cols=colnames(x0)[1:2], alpha=.5)+ # length(unique(mEta$sOurce))
    geom_scatterpie_legend(x0$radius, x=-150, y=-20)+
    scale_y_continuous(limits = c(-60,90), expand = c(0, 0))+
    scale_fill_manual(values=cOl)
ggsave(paste0(pT[2],"DataSrc--map.pdf"), plot = p0, width = 70, height = 30, units = "cm")

x0$radius = x0$radius/5
p0 = ggplot()+theme_bw()+coord_equal()+
    geom_point(data = map.world, aes(x=long, y=lat), cex=.1, color="#00000033")+
    geom_label_repel(aes(x=long, y=lat, label=cOuntry), data=x0, size=7, max.overlaps = getOption("ggrepel.max.overlaps", default = 40))+
    geom_scatterpie(aes(x=long, y=lat, group=cOuntry, r=radius), data=x0, cols=colnames(x0)[1:2], alpha=.5)+
    geom_scatterpie_legend(x0$radius, x=0, y=70)+
    scale_y_continuous(limits = c(25,80), expand = c(0, 0))+
    scale_x_continuous(limits = c(-15,30), expand = c(0, 0))+
    scale_fill_manual(values=cOl)
ggsave(paste0(pT[2],"DataSrc--mapEU.pdf"), plot = p0, width = 30, height = 30, units = "cm")
