#!/bin/env Rscript
# author: ph-u
# script: dNdS_median-cir.r
# desc: plot circle median dN/dS by source
# in: Rscript dNdS_median-cir.r
# out: res/dNdS_median-cir.pdf
# arg: 0
# date: 20240705

cat("Load env:", date(), "\n")
source("p_src.r")
source("p_metabolism_PAO1.r")

#suppressPackageStartupMessages(library(circlize))
library(ggplot2); library(rlang);library(lattice)

d0 = c("Cons", "Neut", "Vari", "No Data")
d0 = data.frame(tYpe = d0, index = c(1:(length(d0)-1),NA), cOl = cBp[1:length(d0)], cex = c(4,12,14,.1))

cat("Import:", date(), "\n")
dNdS.med = read.csv(paste0(pT[1],"dNdS_median-PCA.csv"), header = T, row.names = 1)[,c(2,4,3,5)]
colnames(dNdS.med) = paste(1:ncol(dNdS.med),colnames(dNdS.med), sep = ".")
d.plt = data.frame(sRc = rep(colnames(dNdS.med), each = nrow(dNdS.med)), dNdS = unname(unlist(dNdS.med)))

tHs = median(dNdS.med[dNdS.med>1 & dNdS.med!= Inf], na.rm=T) # .6
tHs = c(tHs, tHs^(-1))^(-1)
dNdS.med0 = dNdS.med[,c(1,3)]
d.plt0 = data.frame(sRc = rep(colnames(dNdS.med0), each = nrow(dNdS.med0)), dNdS = unname(unlist(dNdS.med0)))
tHs0 = median(dNdS.med0[dNdS.med0>1 & dNdS.med0!= Inf], na.rm=T) # .6
tHs0 = c(tHs0, tHs0^(-1))^(-1)
#tHs = c(summary(dNdS.med[dNdS.med>0 & dNdS.med<=1], na.rm=T)[5], median(dNdS.med[dNdS.med>1 & dNdS.med!= Inf], na.rm=T)) # medians between 0-1, 1-max
d.plt[,2] = ifelse(is.na(d.plt[,2]), d0$index[4], ifelse(is.infinite(d.plt[,2]), d0$index[3], ifelse(d.plt[,2] <= tHs[1], d0$index[1], ifelse(d.plt[,2] <= tHs[2], d0$index[2], d0$index[3]))))
d.plt0[,2] = ifelse(is.na(d.plt0[,2]), d0$index[4], ifelse(is.infinite(d.plt0[,2]), d0$index[3], ifelse(d.plt0[,2] <= tHs0[1], d0$index[1], ifelse(d.plt0[,2] <= tHs0[2], d0$index[2], d0$index[3]))))
cat("Selection dN/dS thresholds:\n- conserved: <=",tHs[1],"\n- Neutral/Drift:",tHs[1],"< dN/dS <=", tHs[2], "\n- Variable: >",tHs[2], "\n")
cat("Selection dN/dS thresholds (pwCF & env only):\n- conserved: <=",tHs0[1],"\n- Neutral/Drift:",tHs0[1],"< dN/dS <=", tHs0[2], "\n- Variable: >",tHs0[2], "\n")

for(i in 1:ncol(dNdS.med)){d.plt[which(d.plt[,1] == colnames(dNdS.med)[i]),2] = d.plt[which(d.plt[,1] == colnames(dNdS.med)[i]),2] + .1*(i-2)};rm(i)
for(i in 1:ncol(dNdS.med0)){d.plt0[which(d.plt0[,1] == colnames(dNdS.med0)[i]),2] = d.plt0[which(d.plt0[,1] == colnames(dNdS.med0)[i]),2] + .1*(i-2)};rm(i)
d.tmp = dNdS.med0

d.plt$gEne = row.names(dNdS.med)
d.plt[,2] = ifelse(d.plt$gEne %in% GTF$NCBI.Locus.Tag[which(GTF$Feature.Type=="CDS")], d.plt[,2], NA) # ORF only
d.plt0$gEne = row.names(dNdS.med0)
d.plt0[,2] = ifelse(d.plt0$gEne %in% GTF$NCBI.Locus.Tag[which(GTF$Feature.Type=="CDS")], d.plt0[,2], NA) # ORF only

cat("Plotting:", date(), "\n")
p0 = ggplot() + theme_minimal() + coord_radial(r.axis.inside = T, inner.radius = .3) +
#    xlab("Gene PA number on PAO1") + ylab("Type of dN/dS sequences") +
    xlab("") + ylab("") + ylim(0,3) +
    scale_color_manual(values = set_names(cBp[1:ncol(dNdS.med)], colnames(dNdS.med)), name = "Sampling Source") +
#    scale_x_continuous(breaks = 1:nrow(dNdS.med), label = ifelse(((1:nrow(dNdS.med))%%477)==1, row.names(dNdS.med), "")) +
    scale_x_continuous(breaks = as.numeric(sub("PA","",Mbp$locusTag)), label = Mbp$loc) +
    scale_y_continuous(breaks = 1:max(round(d.plt[,2]), na.rm = T), label = rep("",length(d0$tYpe[-nrow(d0)]))) +
    theme(panel.grid = element_blank(), plot.margin = unit(c(-1,0,-1,0), "cm"), axis.text = element_text(size = 24), axis.text.y = element_text(angle = 45)) +
    geom_point(aes(x = rep(1:nrow(dNdS.med), ncol(dNdS.med)), y = d.plt[,2], group = d.plt[,1], color = d.plt[,1]))
#    geom_point(aes(x = 1:nrow(d.plt), y = 1, shape = d0$shape[match(d.plt[,1], d0$index)], size = d0$cex[match(d.plt[,1], d0$index)]))
ggsave(paste0(pT[2],"dNdS_median-cir.pdf"), plot = p0, width = 10, height = 10)
cat("Plot done:", date(), "\n")
d.plt[,2] = round(d.plt[,2])
write.csv(d.plt[,c(3,1,2)], paste0(pT[1],"dNdS_median-cir.csv"), quote = F, row.names = F)

p0 = ggplot() + theme_minimal() + coord_radial(r.axis.inside = T, inner.radius = .3) +
    xlab("") + ylab("") + ylim(0,3) +
    scale_color_manual(values = set_names(cBp[1:ncol(dNdS.med)], colnames(dNdS.med)), name = "Sampling Source") +
    scale_x_continuous(breaks = as.numeric(sub("PA","",Mbp$locusTag)), label = Mbp$loc) +
    scale_y_continuous(breaks = 1:max(round(d.plt[,2]), na.rm = T), label = rep("",length(d0$tYpe[-nrow(d0)]))) +
    theme(panel.grid = element_blank(), plot.margin = unit(c(-1,0,-1,0), "cm"), axis.text = element_text(size = 24), axis.text.y = element_text(angle = 45)) +
    geom_point(aes(x = rep(1:nrow(dNdS.med0), length(unique(d.plt0[,1]))), y = d.plt0[,2], group = d.plt0[,1], color = d.plt0[,1]))
ggsave(paste0(pT[2],"dNdS_median-cir2.pdf"), plot = p0, width = 10, height = 10)

##### Neutral in CF, non-neutral in env & vice versa #####
dNdS.med0 = dNdS.med[which(row.names(dNdS.med) %in% GTF$Locus.Tag[which(GTF$Feature.Type=="CDS")]),c(1,3)] # ORF, CF vs env only
for(i in 1:ncol(dNdS.med0)){
 dNdS.med0[,i] = ifelse(dNdS.med0[,i] > tHs[2], 1, ifelse(dNdS.med0[,i] < tHs[1], -1, 0))
};rm(i)
dNdS.med0$gNam = GTF$gNam[match(row.names(dNdS.med0),GTF$Locus.Tag)]
dNdS.med0$pRod = GTF$Product.Description[match(row.names(dNdS.med0),GTF$Locus.Tag)]
dNdS.med0$cOnf = GTF$Product.Name.Confidence[match(row.names(dNdS.med0),GTF$Locus.Tag)]

d.med1 = vector(mode = "list", length = 9)
i0=1; for(i in -1:1){for(j in -1:1){
    d.med1[[i0]] = dNdS.med0[which(dNdS.med0[,1]==i & dNdS.med0[,2]==j),]
    i0 = i0+1
}};rm(i,j,i0)

#nrow(dNdS.med0[which(dNdS.med0[,1]<0 & dNdS.med0[,2]<0),])

dNdS.heat = c();for(i0 in -1:1){for(i1 in -1:1){dNdS.heat = c(dNdS.heat,nrow(dNdS.med0[which(dNdS.med0[,1]==i1 & dNdS.med0[,2]==i0),]))}};rm(i0,i1)
dNdS.heat[dNdS.heat>200] = 200
pdf(paste0(pT[2],"dNdS_median-heat.pdf"), width = 8.27, height = 11.69, paper = "a4")
print(levelplot(matrix(dNdS.heat, nrow = 3), scales=list(x=list(rot=90)), col.regions = colorRampPalette(c("blue", "red"))(100), xlab = "Environmental", ylab = "pwCF"))
invisible(dev.off())

# i=3;i = d.med1[[i]][which(d.med1[[i]]$cOnf=="Class 1"),];dim(i);paste0(row.names(i), collapse = ",");paste0(paste0(row.names(i)," (",i$gNam,")"), collapse = ", ");rm(i)
##### f: reformat PAnum and gene names #####
gNamReformat = function(i){
    i0 = dNdS.med0[which(row.names(dNdS.med0) %in% strsplit(i,"\n")[[1]]),]
    i1 = i0[which(i0$cOnf=="Class 1"),]
    return(paste0(paste0(row.names(i1)," [",i1$gNam,"]"), collapse = ", "))}

##### Get differential-selected ORFs #####
d.diffORFs = dNdS.med0
d.diffORFs$pRod = gsub(",",";",d.diffORFs$pRod)
d.diffORFs$dNdS.median.env = d.diffORFs$dNdS.median.CF = NA
src.diffORFs = list(mEta$assemblyInfo.genbankAssmAccession[mEta$sOurce=="Cystic fibrosis"], mEta$assemblyInfo.genbankAssmAccession[mEta$sOurce=="Environmental"])
for(i in 1:nrow(d.diffORFs)){
  i0 = list.files(pT[3],paste0(row.names(d.diffORFs)[i],"--"), full.names=T)
  i0 = read.csv(i1 <- i0[grep("dbSum",i0)], header = T)
  cat(i1,"     \r")
  i0$clinical = read.table(text = sub("_ASM","@",i0$clinical), sep = "@")[,1]
  i0$dNdS = abs(i0$dNdS)
  d.diffORFs[i,ncol(d.diffORFs)+-1:0] = c(median(i0$dNdS[which(i0$clinical %in% src.diffORFs[[1]])], na.rm = T), median(i0$dNdS[which(i0$clinical %in% src.diffORFs[[2]])], na.rm = T))
};rm(i,i0,i1);cat("\n")
write.csv(d.diffORFs, paste0(pT[2],"differentialORF--all.csv"), row.names = T, quote = F)
d.diffORFs = d.diffORFs[which(!is.na(d.diffORFs[,1]) & !is.na(d.diffORFs[,2]) & !(d.diffORFs[,1]<0 & d.diffORFs[,2]<0)),-ncol(d.diffORFs)]
write.csv(d.diffORFs, paste0(pT[2],"differentialORF.csv"), row.names = T, quote = F)
