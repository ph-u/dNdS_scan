#!/bin/env Rscript
# author: ph-u
# script: fig_vgrG5.r
# desc: plot VgrG5 operon
# in: Rscript fig_vgrG5.r
# out: res/fig_vgrG5.jpeg
# arg: 0
# date: 20250602

source("p_src.r")

gDF = GTF[grep("PA5086",GTF$Locus.Tag):grep("PA5090",GTF$Locus.Tag),c("Coordinates","start","end")]

jpeg(paste0(pT[2],"fig_vgrG5.jpeg"), width = 1500, height = 1000, res = 300)
par(mar = c(0,0,0,0))
plot(x = c(gDF$start,gDF$end), y = rep(0,nrow(gDF)*2), xlim = range(c(gDF$start,gDF$end)), ylim = c(-1,1), col = "#00000000", xlab = "", ylab = "")
arrows(x0 = gDF$end, x1 = gDF$start, y0 = rep(0,nrow(gDF)), y1 = rep(0,nrow(gDF)), length = .1, angle = 20)
text(x = rowMeans(gDF[,-1]), y = .1, labels = paste0("PA50",86:90), cex = .8)
invisible(dev.off())
