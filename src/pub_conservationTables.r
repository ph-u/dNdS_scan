#!/bin/env Rscript
# author: ph-u
# script: pub_conservationTables.r
# desc: get ORF conservation tables for categorizations and volcano plot
# in: Rscript pub_conservationTables.r
# out: res/p_conservationTables.csv
# arg: 0
# date: 20250322

##### set env #####
source("p_src.r")
library(EnhancedVolcano)
diffORF = read.csv(paste0(pT[2],"differentialORF--all.csv"), row.names = 1, header = T)
diffORF$gene = row.names(diffORF)
vOlcano = read.csv(paste0(pT[1],"dNdS_infect-vol.csv"), header = T)
vOlcano = vOlcano[which(vOlcano$cond1=="Cystic fibrosis" & vOlcano$cond2=="Environmental" & vOlcano$gene %in% GTF$Locus.Tag[GTF$Feature.Type=="CDS"]),]

##### Get Volcano plot color details #####
aaa = EnhancedVolcano(vOlcano,
                lab = vOlcano$gene,
                x = "log2FC",
                y = 'p.val',
                xlab = bquote(~Log[2] ~ "[" ~ Log[10] ~ "median ratio] change"),
                title = paste0(capFirst(unique(vOlcano$cond1))," (L) vs ",capFirst(unique(vOlcano$cond2))," (R)"),
                pCutoff= .1, #FDR 
                pCutoffCol = "p.adj", # EnhancedVolcano suggestion: false-discovery rate
                FCcutoff = 1,
                pointSize = 2,
                labSize = 0,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendLabels = c('NS', expression(Log[2] ~ LMRC), 'Adjusted p-value', expression(adjusted ~ p - value ~ and ~ log[2] ~ LMRC)), 
                legendPosition = 'bottom',
                legendLabSize = 12, 
                legendIconSize = 4,
                drawConnectors = F,
                widthConnectors = 0)

vOlcano$log10P = -log10(vOlcano$p.val)
aaa$data$Sig = as.character(aaa$data$Sig)
vOlcano$volcano.cOlor = aaa$data$Sig[match(vOlcano$gene,aaa$data$gene)]

d0 = merge(vOlcano,diffORF, by="gene", all=T)
d0 = d0[,c(1,11,14,15,9,10,4:8,12)]
rm(diffORF, vOlcano, aaa)

##### Get Count Tables #####
d0$X1.Cystic.fibrosis[is.na(d0$X1.Cystic.fibrosis)] = "NA"
d0$X3.Environmental[is.na(d0$X3.Environmental)] = "NA"
d0$volcano.cOlor = ifelse(d0$volcano.cOlor=="FC_P","Red",ifelse(d0$volcano.cOlor=="FC","Green",ifelse(d0$volcano.cOlor=="P","Blue","Grey")))
print(table(d0[,c("X1.Cystic.fibrosis","X3.Environmental")]))
print(table(d0$volcano.cOlor))

d0 = d0[,c("gene","gNam","dNdS.median.CF","dNdS.median.env","X1.Cystic.fibrosis","X3.Environmental","log2FC","p.val","p.adj","log10P","volcano.cOlor","pRod")]
write.csv(d0, paste0(pT[1],"p_conservationTables.csv"), row.names = F, quote = F)
# a0 = a[which((a$log2FC==Inf) & a$log10P>-log10(.1)),];a0
