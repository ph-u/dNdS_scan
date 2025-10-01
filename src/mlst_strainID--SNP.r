#!/bin/env Rscript
# author: ph-u
# script: mlst_strainID--SNP.r
# desc: Calculate SNP count for isolates suspecting to be the same strain
# in: Rscript mlst_strainID--SNP.r
# out: res/mlst_strainID--SNP.csv
# arg: 0
# date: 20250423

source("p_src.r")
set.seed(8389)
suppressMessages(library(igraph))
mlst.st = read.csv(paste0(pT[2],"mlst_REALPHY--metadata.csv"), header = T)
mlst.id = read.csv(paste0(pT[2],"mlst_strainID.csv"), header = T, row.names = 1)
gLen = read.csv(paste0(pT[1],"mlst_strainID--gLen.csv"), header = T)
st.uniq = unique(mlst.st$seqType) # c(155,179,27)

##### Comparison pairwise set #####
i4 = c();for(i in 1:length(st.uniq)){
  i1 = mlst.st[which(mlst.st$seqType==st.uniq[i]),1]
  for(i2 in 1:(length(i1)-1)){for(i3 in (i2+1):length(i1)){i4 = c(i4,paste0(i1[i2]," vs ",i1[i3]))}}
};rm(i)

##### Calculate total SNP differences
i0 = c();for(i in 1:length(st.uniq)){
  i0 = c(i0,colSums(read.csv(paste0(pT[1], "mlst_stSNP--ST",st.uniq[i],".csv"), header = T, row.names = 1), na.rm = T))
};rm(i)

write.csv(data.frame(group = colnames(mlst.id), SNP = i0, sameStrain = ifelse(i0<25e3,"Y","")), paste0(pT[2],"mlst_strainID--SNP.csv"), row.names = F, quote = F)

##### Build network graph #####
net0 = data.frame(group = colnames(mlst.id), pair = i4, SNP = i0, sim = 1-i0/sum(gLen$gLen))
net0$scaled = ifelse(net0$SNP>25e3,0,ifelse(net0$SNP<15e3,0,1))

netNam = data.frame(UserID = 1:nrow(mlst.st), Name = mlst.st$assemblyInfo.genbankAssmAccession, ST = mlst.st$seqType, src = mlst.st$sOurce)
net0 = cbind(net0,read.table(text=gsub(" vs ","!",net0$pair), sep = "!"))
net0$V1 = netNam$UserID[match(net0$V1,netNam$Name)]
net0$V2 = netNam$UserID[match(net0$V2,netNam$Name)]

gDF = graph_from_data_frame(net0[,c("V1","V2")], directed = F, vertices = netNam)
V(gDF)$color = ifelse(V(gDF)$ST=="155",cBp[2],ifelse(V(gDF)$ST=="27",cBp[6],cBp[4]))
E(gDF)$color = ifelse(net0$SNP<25e3,"#000000ff", "#00000000") # black = clonal/strain
E(gDF)$lty = ifelse(net0$SNP>15e3,1,2) # dash = clonal
V(gDF)$shape = ifelse(V(gDF)$src=="Cystic fibrosis", "circle","square")

pdf(paste0(pT[2],"mlst_strainID--net.pdf"), width = 10, height = 10)
# https://r-graph-gallery.com/248-igraph-plotting-parameters.html
plot(gDF, layout = layout_with_fr(gDF), vertex.label = V(gDF)$Name, vertex.shape = V(gDF)$shape,
  vertex.frame.color = "white", vertex.size = 10, vertex.label.family="Helvetica",
  edge.color = E(gDF)$color, edge.lty = E(gDF)$lty)
invisible(dev.off())
