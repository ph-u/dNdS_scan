#!/bin/env Rscript
# author: ph-u
# script: varTypeBar_circle.r
# desc: sequence variation type in all PAO1 dN/dS processed genes
# in: Rscript varTypeBar_circle.r
# out: res/vt_Bar_c.pdf
# arg: 0
# date: 20240320, 20240426

# https://r-graph-gallery.com/295-basic-circular-barplot.html
# https://r-graph-gallery.com/circular-barplot.html
# https://r-graph-gallery.com/web-circular-lollipop-plot-with-ggplot2.html
# https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html
source("p_src.r")
library(gridExtra); library(rlang)

fTag = unique(mEta$sOurce)[order(unique(mEta$sOurce))]
f.count = table(mEta$sOurce)
f.title = paste(names(f.count)," n =",f.count, sep = " ")
vtBar = read.csv(paste0(pT[2],"vt_Bar.csv"), header = T, row.names = 1)

##### Set record format #####
f.r0 = vector(mode = "list", length = length(fTag))
names(f.r0) = fTag
for(i in 1:length(f.r0)){
    i0 = as.data.frame(matrix(0, nr = nrow(f), nc = ncol(vtBar)))
    row.names(i0) = f$gEne
    colnames(i0) = colnames(vtBar)
    f.r0[[i]] = i0
};rm(i, i0)

cat("Categorizing variation types:",date(),"\n")
for(i0 in 1:nrow(f)){
    cat("Processing",f$gEne[i0],":",date(),"     \r")
    eSx = unique(read.csv(paste0(pT[3],paste0("PAO1_107_",f$gEne[i0],"--dbSum.csv")), header = T))
    eSx$sRc = mEta$sOurce[match(read.table(text = gsub("_ASM","@",eSx$clinical), sep = "@")[,1], mEta$assemblyInfo.genbankAssmAccession)]

    vType = strsplit(eSx$varType, ":")
    eSx$vType = eSx$varType
    for(i in 1:length(vType)){if(length(vType[[i]]) > 1){eSx$vType[i] = paste(vType[[i]][c(1,length(vType[[i]]))], collapse = ".")}};rm(i, vType)

    eSx.0 = as.data.frame.matrix(table(eSx[,c("sRc","vType")]))/f.count*100
    for(i in 1:nrow(eSx.0)){f.r0[[i]][i0,match(colnames(eSx.0), colnames(f.r0[[1]]))] = eSx.0[i,]};rm(i)
};rm(i0,eSx, eSx.0)
for(i in 1:length(f.r0)){f.r0[[i]][is.na(f.r0[[i]])] = 0};rm(i)
save(f.r0, file=paste0(pT[1],"varTypeBar_circle.rda"))

cat("\nPlotting:",date(),"\n")
xLab = Mbp$loc[match(row.names(f.r0[[1]]),Mbp$locusTag)]
xLab[is.na(xLab)] = ""
# https://stackoverflow.com/questions/65137163/geom-bar-removed-rows-containing-missing-values-but-doesnt
# https://www.tidyverse.org/blog/2024/03/ggplot2-3-5-0-coord-radial/
for(i in 1:length(fTag)){
    pCir = ggplot() +
        geom_col(aes(x = as.factor(rep(row.names(f.r0[[i]]), ncol(f.r0[[i]]))), y = unname(unlist(f.r0[[i]])), fill = rep(colnames(f.r0[[i]]), each = nrow(f.r0[[i]]))), alpha=1) +
        scale_fill_manual(values = set_names(cBp[1:ncol(f.r0[[1]])], colnames(f.r0[[1]])), name = "Type") +
        ggtitle(f.title[i]) + xlab("") + ylab("") + #xlab("Gene PA number on PAO1") + ylab("Proportion of the Type of Sequence Variation (%)") +
#        scale_x_discrete(label = ifelse(((1:nrow(f.r0[[i]]))%%477)==1, row.names(f.r0[[i]]), "")) +
        scale_x_discrete(label = xLab) +
        theme(axis.text = element_text(size = 24),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.margin = unit(rep(-1,4), "cm")) +
        theme_minimal() + coord_radial(r.axis.inside = T, inner.radius = .1) # ggplot2 >=3.5.0 (3.5.1)
    ggsave(paste0(pT[2],"vt_Bar_c_",gsub(" ","-",fTag[i]),".pdf"), plot = pCir, width = 7, height = 7)
    cat("Exported barplot",fTag[i],":",date(),"\n")
};rm(i)
rm(xLab)

##### Contrast type of variation between CF and env #####
#cat("Contrast varType: CF vs env -- ",date(),"\n")

#vtContrast = f.r0[[1]]/f.r0[[2]]
#vtContrast[is.na(vtContrast)] = 1
#for(i in 1:ncol(vtContrast)){vtContrast[is.infinite(vtContrast[,i]),i] = max(unlist(vtContrast)[is.finite(unlist(vtContrast))])};rm(i)

#jpeg(paste0(pT[2],"vt_Bar_noHit.jpeg"), width = 2400, height = 1200, res = 300)
#par(mar = c(2,0,0,0)+.1)
#plot(x = 1:nrow(vtContrast), y = log10(vtContrast[,5]), type = "p", pch = 3, xaxt = "n", yaxt = "n", ylab = "", xlab = "", cex = .2)
#axis(1, at = as.numeric(sub("PA","",Mbp$locusTag)), labels = Mbp$loc)
#abline(v = c(2125,2384), col = "#ff00ffff")
#invisible(dev.off())

#cat("Contrast done -- ",date(),"\n")

#pdf(paste0(pT[2],"vt_Bar_circle.pdf"))
#grid.arrange(pCir[[1]],pCir[[2]],pCir[[3]],pCir[[4]],pCir[[5]], ncol = 1)
#invisible(dev.off())
# https://stackoverflow.com/questions/1249548/side-by-side-plots-with-ggplot2#3935554
