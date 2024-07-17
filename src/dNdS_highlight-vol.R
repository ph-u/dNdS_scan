#!/bin/env Rscript
# author: ph-u
# script: dNdS_highlight-vol.R
# desc: volcano plot gene highlight
# in: Rscript dNdS_highlight-vol.R
# out: dNdS_highlight-vol.pdf
# arg: 0
# date: 20240714

#pOri = getwd()
#setwd("../../../1_02_wholeGenomeDNDS/src")
source("p_src.r")
vOl = read.csv(paste0(pT[1],"dNdS_infect-vol.csv"), header = T)
eSs = read.table(paste0(pT[4],"pnas.1419677112.sd01.csv"), sep = "\t", header = T)
eSs = eSs[,c(1,grep("Essen", colnames(eSs)))]
eSs[eSs=="Essential"] = 0
eSs[eSs=="Nonessential"] = 1
for(i in 2:ncol(eSs)){eSs[,i] = as.numeric(eSs[,i])};rm(i)
#setwd(pOri);rm(pOri)

hL = paste0("PA",c("0424","0425","0762","0763","0958","0861","1288","1430","1798",2020,2426,2491,2492,3168,3458,3477,3545,3974,4266,4367,4522,4601,5291), sep = "") # c("0674", 2282:2286) # c("0674", 4984, 5402:5407)

library(EnhancedVolcano)

#d.t = vOl[grep("Cystic", vOl$cond1),]
d.t = vOl[which(vOl$cond1=="Cystic fibrosis" & vOl$cond2=="Environmental"),]
#print(d.t$gene[is.na(d.t$log2FC)])
d.t = d.t[!is.na(d.t$log2FC),]
d.t$log2FC[is.infinite(d.t$log2FC)] = ifelse(d.t$log2FC[is.infinite(d.t$log2FC)]>0,max(d.t$log2FC[is.finite(d.t$log2FC)])+1, min(d.t$log2FC[is.finite(d.t$log2FC)])-1)

pdf(paste0(pT[2],"dNdS_volHigh_",hL[1],"_",unique(d.t$cond1),"_",unique(d.t$cond2),"-vol.pdf"), width = 10, height = 12)
print(EnhancedVolcano(d.t,
    lab = ifelse(d.t$gene %in% hL, d.t$gene,""),
    x = "log2FC",
    y = 'p.val',
    xlab = bquote(~Log[2] ~ "[" ~ Log[10] ~ "median ratio] change"),
    title = paste0(capFirst(unique(d.t$cond1))," (L) vs ",capFirst(unique(d.t$cond2))," (R)"),
    pCutoff= .1, #FDR 
    pCutoffCol = "p.adj", # EnhancedVolcano suggestion: false-discovery rate
    FCcutoff = 1,
    pointSize = 2,
    labSize = 4,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    legendLabels = c('NS', expression(Log[2] ~ LMRC), 'Adjusted p-value', expression(adjusted ~ p - value ~ and ~ log[2] ~ LMRC)),
    legendPosition = 'bottom',
    legendLabSize = 12,
    legendIconSize = 4,
    drawConnectors = T,
    max.overlaps = Inf,
    widthConnectors = .5))
invisible(dev.off())

d.t$p.LOG10 = -log10(d.t$p.val)
pdf(paste0(pT[2],"dNdS_volHigh_",hL[1],"_",unique(d.t$cond1),"_",unique(d.t$cond2),"-man.pdf"), width = 50, height = 10)
par(mar = c(5,6,1,1)+.1, cex.lab = 3, cex.axis = 3)
plot(x = 1:nrow(d.t), y = d.t$p.LOG10,
# x = extreme +FC (variable in CF), + = extreme -FC (conserved in CF), triangle = dN/dS indifferent between CF & env
    pch = ifelse(d.t$log2FC > 1, 4, ifelse(d.t$log2FC < (-1), 3, 2)),
    xlab = "Genomic position", ylab = bquote(~-Log[10] ~ italic(P)),
# big = essential genes (https://doi.org/10.1073/pnas.1419677112)
    cex = ifelse(d.t$gene %in% eSs$Locus.ID[rowSums(eSs[,-1])==0], 5,2),
# orange = high mutational burden in CF (https://doi.org/10.1126/science.adi0908)
    col = ifelse(d.t$gene %in% hL, cBp[3],paste0(substr(cBp[1],1,7),"55")))
# pink text = gene with significance beyond threshold
text(x = 1:nrow(d.t), y = d.t$p.LOG10+.3, labels = ifelse(d.t$p.LOG10>=9, d.t$gene,""), cex = 1.2, col = cBp[9])
abline(h = -log10(.1), lty = 2, lwd = 2, col = cBp[9])
invisible(dev.off())
write.csv(d.t, paste0(pT[1],"dNdS_volHigh_",hL[1],"_",unique(d.t$cond1),"_",unique(d.t$cond2),"-man.csv"), row.names = F, quote = F)
