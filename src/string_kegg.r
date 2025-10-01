#!/bin/env Rscript
# author: ph-u
# script: string_kegg.r
# desc: Lint KEGG biopathway map
# in: Rscript string_kegg.r
# out: NA
# arg: 0
# date: 20250501

source("p_src.r")
library(stringr)
#library(pathview) # https://rdrr.io/bioc/pathview/man/download.kegg.html

cLus = read.csv(paste0(pT[2],"string_diff.csv"), header = T)
cLus$gP = paste(cLus$cID,cLus$CFenv, sep = "_")

cPath = as.data.frame(matrix(nrow = length(unique(cLus$gP)), ncol = 3))
colnames(cPath) = c("cluster", "pathway", "kegg")
cPath$cluster = unique(cLus$gP)
for(i in 1:nrow(cPath)){
  i0 = pAthway[which(pAthway$Locus.Tag %in% cLus$PAnum[which(cLus$gP==cPath$cluster[i])]),c("Locus.Tag","Pathway.Name","Pathway.Xref","Pathway.Database")]
  i1 = iNterpro[which(iNterpro$Locus.Tag %in% cLus$PAnum[which(cLus$gP==cPath$cluster[i])]),c("Locus.Tag","Product","External.DB")]
  cPath[i,-1] = c(gsub(",","",paste0(c(i0$Pathway.Name,i1$Product), collapse = ";")),sub("^;","",gsub(",","",paste0(unique(i0$Pathway.Xref)[order(unique(i0$Pathway.Xref))], collapse = ";"))))
};rm(i,i0,i1)
write.csv(cPath, paste0(pT[1],"string_kegg.csv"), row.names = F, quote = F)

##### Keywords aggregation #####
cPath$pathway = gsub("< sup>","",gsub("<sup>","",gsub("[;():./]"," ",tolower(cPath$pathway))))

for(i in c(
"probable","pseudomonas","putative","unsaturated","aeruginosa","replicative","superfamily","superpathway","conserved","and"
)){ cPath$pathway = gsub(i,"",cPath$pathway) }

for(i in c(
"cell","aromatic","50s","rna","aerobic","dna","weight","general","osmotically","efflux","immunity","secretion","chain","secondary","response","to","i","abc","common","diverse","electron","fad","major","metaboli","multidrug","nad","nucleotide","organic","other","al","d1","de","novo","backbone"
)){ cPath$pathway = gsub(paste0(i," "),paste0(i,"-"),cPath$pathway) }

for(i in c(
"s5","carrier","coa","salvage","blocks","unit","synthesis","subunit","coupling","kinase","fixation","excision","vi","genes","acid","1a","factor","membrane","cyanide","translocation","system","signal","sensing","sensor","shock","secretion","utilization","response","reaction","repair","resistance","to","formation","assimilation","binding","i","degradation","export","hii","ligase","mutase","oxi","path","recomb","trans","fermentation","cycle","coenzyme","cat","bypass"
)){ cPath$pathway = gsub(paste0(" ",i),paste0("-",i),cPath$pathway) }

for(i in c(
"type","abc","ferredoxin","l-asparaginase","beta-lactam","bacterial-secretion","transcrip","hii-ribonuclease","topoisomerase","oxidoreductase"
)){ cPath$pathway = gsub(paste0("-",i),paste0(" ",i),cPath$pathway) }

cPath$pathway = gsub("regulator-","regulator ",gsub("biosynthesis-","biosynthesis ",sub("-$","",cPath$pathway)))

rm(i)
pRaw = unique(unlist(strsplit(cPath$pathway, " ")))
pRaw = pRaw[nchar(pRaw)>1]
pRaw = pRaw[-which(nchar(pRaw)==4)]
pRaw = pRaw[order(pRaw)]
pRaw = pRaw[-which(pRaw %in% c("--","-acidic","-formation","-path","[positive","a-pathway","acetyl-coa-assimilation","activation","aer","amino-acid","anaerobic-respiration","apr","azor3","bett3","biotin-synthesis","by","chain","chemotaxis]","class-iii-","component","diverse-environments","ferredoxin-i-","glycolysis-i-","glycolysis-ii-","hii-ribonuclease-hii-","i-","ii-","l6","metabolic-pathways","metabolism","microbial-metabolism-in","of","one","phosphorylation-oxidative","propionate-catabolic","protease-secretion-protein","protein","protein-1a","protein-in","reductase-subunit","riboflavin-kinase","rpmj2","spee2","tagf1","tagt1","tca-cycle-i-","thiol","tli5b3","tsse1","tssj1","tca-cycle","synthase-iii-3-oxoacyl-[acyl-carrier-protein]","synthase-iii-nosf","serine-isocitrate","for"))]

i0 = c();for(i in 1:length(pRaw)){if(length(grep(paste0(pRaw[i],"s"),pRaw, fixed = T))>0){i0 = c(i0,grep(paste0(pRaw[i],"s"),pRaw, fixed = T))}};rm(i)
pRaw = pRaw[-unique(i0)];rm(i0)

##### PCA #####
## cluster & CF/env
pcaRaw = as.data.frame(matrix(nrow = nrow(cPath), ncol = length(pRaw)+3))
pcaRaw[,1] = cPath$cluster
pcaRaw[,2:3] = matrix(unlist(strsplit(cPath$cluster, "_")), ncol = 2, byrow = T)

for(i in 1:length(pRaw)){
  pcaRaw[,i+3] = str_count(cPath$pathway,pattern = fixed(pRaw[i]))
};rm(i)

pCa = prcomp(pcaRaw[,-(1:3)], center = T, scale. = T)

i=3 #for(i in 2:3){
  g0 = ggbiplot(pCa, var.scale = 1, labels = pcaRaw[,1], groups = pcaRaw[,i], ellipse = T, ellipse.prob = .95, labels.size = 20, var.axes = F) + theme_bw()+
    coord_cartesian(xlim = c(-3, 6)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_color_manual(values=setNames(cBp[1:length(unique(pcaRaw[,i]))], unique(pcaRaw[,i])))
  ggsave(paste0(pT[2],"string_kegg--pca",ifelse(i==2,"Clus","CFev"),".jpeg"), plot = g0, width = 10, height = 10)
#};rm(i)

## cluster only
pcaRaw[nrow(pcaRaw)+1,] = 0
pcaRaw0 = pcaRaw[seq(1,nrow(pcaRaw),2),-c(1:3)] + pcaRaw[seq(2,nrow(pcaRaw),2),-c(1:3)]
row.names(pcaRaw0) = 1:nrow(pcaRaw0)
pcaRaw = pcaRaw[-nrow(pcaRaw),]

pCa = prcomp(pcaRaw0, center = T, scale. = T)

g0 = ggbiplot(pCa, var.scale = 1, labels = paste0("G",row.names(pcaRaw0)), labels.size = 20, var.axes = F) + theme_bw()+
  coord_cartesian(xlim = c(-5, 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste0(pT[2],"string_kegg--pcaClus.jpeg"), plot = g0, width = 10, height = 10)

##### export #####
pcaRaw0 = as.data.frame(t(pcaRaw))
pcaRaw0$text = c(rep("",3),pRaw)
write.csv(pcaRaw0[,c(ncol(pcaRaw0),1:(ncol(pcaRaw0)-1))], paste0(pT[1],"string_kegg--pcaRaw.csv"), row.names = F, quote = F)

##### KEGG pathway plots data #####
# https://www.genome.jp/kegg/mapper/color.html
cBpP = cBp[-c(1,9)]
cLus$col1 = paste0(substr(cBpP[cLus$cID],1,7),ifelse(cLus$CFenv=="CF","FF","77"))
cLus$col2 = ifelse(cLus$CFenv=="CF","#00ff00ff","#0000ffff")
write.table(cLus[,c("PAnum","col2","col1")], paste0(pT[1],"string_kegg--coPub.csv"), sep = "\t", col.names = F, row.names = F, quote = F)
cLus$col2 = ifelse(cLus$CFenv=="CF","#ff0000ff","#0000ffff")
write.table(cLus[,c("PAnum","col2","col1")], paste0(pT[1],"string_kegg--color.csv"), sep = "\t", col.names = F, row.names = F, quote = F)
write.table(cLus[which(cLus$cID %in% c(6,7,10)),c("PAnum","col2","col1")], paste0(pT[1],"string_kegg--colSel.csv"), sep = "\t", col.names = F, row.names = F, quote = F)
write.table(cLus[-which(cLus$cID %in% c(6,7,10)),c("PAnum","col2","col1")], paste0(pT[1],"string_kegg--colOth.csv"), sep = "\t", col.names = F, row.names = F, quote = F)

## Color Legend
jpeg(paste0(pT[2],"string_kegg--color.jpeg"), width = 1000, height = 1000, res = 300)
par(mar=c(0,0,0,0))
plot(x = c(-3.5,5), y = c(-2.5,2), col = "#00000000")
text(x = rep((1:8)-4, 5), y = rep(2:-2, each = 8)-rep(c(0,.5),20), labels = c(unique(cLus$gP),""), col = c(unique(cLus$col1),"#00000000"), cex = 1)
invisible(dev.off())
