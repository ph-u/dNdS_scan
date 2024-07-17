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
library(ggplot2); library(rlang)

d0 = c("Cons", "Neut", "Vari", "No Data")
d0 = data.frame(tYpe = d0, index = c(1:(length(d0)-1),NA), cOl = cBp[1:length(d0)], cex = c(4,12,14,.1))

cat("Import:", date(), "\n")
dNdS.med = read.csv(paste0(pT[1],"dNdS_median-PCA.csv"), header = T, row.names = 1)[,c(2,5,3,4)]
colnames(dNdS.med) = paste(1:ncol(dNdS.med),colnames(dNdS.med), sep = ".")
d.plt = data.frame(sRc = rep(colnames(dNdS.med), each = nrow(dNdS.med)), dNdS = unname(unlist(dNdS.med)))

tHs = median(dNdS.med[dNdS.med>1], na.rm=T) # .6
tHs = c(tHs, tHs^(-1))^(-1)
#tHs = c(summary(dNdS.med[dNdS.med>0 & dNdS.med<=1], na.rm=T)[5], median(dNdS.med[dNdS.med>1], na.rm=T)) # medians between 0-1, 1-max
d.plt[,2] = ifelse(is.na(d.plt[,2]), d0$index[4], ifelse(is.infinite(d.plt[,2]), d0$index[3], ifelse(d.plt[,2] <= tHs[1], d0$index[1], ifelse(d.plt[,2] <= tHs[2], d0$index[2], d0$index[3]))))
cat("Selection dN/dS thresholds:\n- conserved: <=",tHs[1],"\n- Neutral/Drift:",tHs[1],"< dN/dS <=", tHs[2], "\n- Variable: >",tHs[2], "\n")

for(i in 1:ncol(dNdS.med)){d.plt[which(d.plt[,1] == colnames(dNdS.med)[i]),2] = d.plt[which(d.plt[,1] == colnames(dNdS.med)[i]),2] + .1*(i-2)};rm(i)

cat("Plotting:", date(), "\n")
p0 = ggplot() + theme_bw() + coord_radial(r.axis.inside = T, inner.radius = .3) +
    xlab("Gene PA number on PAO1") + ylab("Type of dN/dS sequences") + ylim(0,3) +
    scale_color_manual(values = set_names(cBp[1:ncol(dNdS.med)], colnames(dNdS.med)), name = "Sampling Source") +
    scale_x_continuous(breaks = 1:nrow(dNdS.med), label = ifelse(((1:nrow(dNdS.med))%%477)==1, row.names(dNdS.med), "")) +
    scale_y_continuous(breaks = 1:max(round(d.plt[,2]), na.rm = T), label = d0$tYpe[-nrow(d0)]) +
    theme(panel.grid = element_blank(), plot.margin = unit(c(-1,0,-1,0), "cm"), axis.text.y = element_text(size = 12)) +
    geom_point(aes(x = rep(1:nrow(dNdS.med), ncol(dNdS.med)), y = d.plt[,2], group = d.plt[,1], color = d.plt[,1]))
#    geom_point(aes(x = 1:nrow(d.plt), y = 1, shape = d0$shape[match(d.plt[,1], d0$index)], size = d0$cex[match(d.plt[,1], d0$index)]))
ggsave(paste0(pT[2],"dNdS_median-cir.pdf"), plot = p0, width = 10, height = 10)
cat("Plot done:", date(), "\n")
d.plt$gEne = row.names(dNdS.med)
d.plt[,2] = round(d.plt[,2])
write.csv(d.plt[,c(3,1,2)], paste0(pT[1],"dNdS_median-cir.csv"), quote = F, row.names = F)
