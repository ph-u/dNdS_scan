#!/bin/env Rscript
# author: ph-u
# script: vt_corTest.r
# desc: categorical sequence variation type correlation test
# in: Rscript vt_corTest.r
# out: NA
# arg: 0
# date: 20240425

# https://towardsdatascience.com/the-search-for-categorical-correlation-a1cf7f1888c9
# https://www.statology.org/correlation-between-categorical-variables/
# https://search.r-project.org/CRAN/refmans/confintr/html/cramersv.html

source("p_src.r")
library(confintr) # v1.0.2
g0 = "PA2185"

##### Seq variation type correlation: Cramer's V #####
dbSum.list = list.files(pT[3],"dbSum") # compare with PA2185
cat("Start Cramer's V calculation:",date(),"\n")
vt_cData = vector(mode = "list", length = length(dbSum.list))
xRf = read.csv(paste0(pT[3],dbSum.list[grep(g0,dbSum.list)]), header = T)
x1 = read.table(text = xRf$varType[grep(":",xRf$varType)], sep = ":")
xRf$varType[grep(":",xRf$varType)] = paste(x1[,1],x1[,ncol(x1)], sep = ":");rm(x1)

vt.cv = as.data.frame(matrix(nr = length(dbSum.list), nc = 1+length(unique(mEta$sOurce))))
colnames(vt.cv) = c("gene", unique(mEta$sOurce)[order(unique(mEta$sOurce))])
vt.cv$gene = read.table(text = gsub("_PA","@PA",gsub("--","@",dbSum.list)), sep = "@")[,2]

for(i in 1:nrow(vt.cv)){
    cat("Processing",vt.cv$gene[i],"-",date(),"     \r")
    x0 = read.csv(paste0(pT[3],dbSum.list[i]), header = T)
    if(length(grep(":",x0$varType))>0){
        x1 = read.table(text = x0$varType[grep(":",x0$varType)], sep = ":")
        x0$varType[grep(":",x0$varType)] = paste(x1[,1],x1[,ncol(x1)], sep = ":")
        rm(x1)
    }
    vt.x = merge(xRf[,c("clinical","varType")], x0[,c("clinical","varType")], by = "clinical", all = T)
    vt.x$sRc = mEta$sOurce[match(read.table(text = gsub("_ASM","@",vt.x$clinical), sep = "@")[,1], mEta$assemblyInfo.genbankAssmAccession)]
    for(i0 in 2:ncol(vt.cv)){
        vt.cv[i,i0] = cramersv(vt.x[which(vt.x$sRc==colnames(vt.cv)[i0]),2:3])
    };rm(i0)
};rm(i, x0);cat("\nCramer's V calculation done:",date(),"\n")
write.csv(vt.cv, paste0(pT[1],"vt_cramersv.csv"), row.names = F, quote = F)

##### Plot genomic correlation distribution by source #####
tHs = .7; cBl = cBp[1:ncol(vt.cv)][-5]
pdf(paste0(pT[2],"vt_cramersv.pdf"), width = 50, height = 21)
par(mar = c(5,5,0,0)+.1, mfrow = c(ncol(vt.cv)-1,1), cex.axis = 2, cex.lab = 2)
for(i in 2:ncol(vt.cv)){
    plot(row.names(vt.cv), vt.cv[,i], col = gsub("FF$","44",cBl[i-1]), type = "l", lty = 1, lwd = 2, xlab = "gene position", ylab = "Cramer's V correlation coefficient")
    abline(h=tHs, lty = 2, col = "#00000077")
    vt.cv0 = vt.cv[vt.cv[,i]>tHs,]
    text(x = row.names(vt.cv0), y = rep(c(.9,.87,.73,.95,.76,.84,.97,.82,.78,.93,.86,.71,1), 100)[1:nrow(vt.cv0)], labels = vt.cv0[,1], col = cBl[i-1])
    legend("topleft", legend = colnames(vt.cv)[-1], pch = 20, cex = 1.2, bg = "#ffffffff", col = gsub("FF$","44",cBl))
};rm(i, vt.cv0)
invisible(dev.off())

vt.cv0 = data.frame(sRc = rep(colnames(vt.cv)[-1], each = nrow(vt.cv)-1), cramersv = unname(unlist(vt.cv[which(vt.cv$gene!=g0),-1])))
pdf(paste0(pT[2],"vt_cramersv_box.pdf"), width = 10, height = 7)
par(mar = c(3,4,0,0)+.1, cex.axis = 1.2)
boxplot(vt.cv0$cramersv ~ as.factor(vt.cv0$sRc), col = "#00000000", pch = 19, xlab = "", ylab = paste0("Cramer's V correlation coefficient against ",g0))
invisible(dev.off())

print(pairwise.wilcox.test(vt.cv0$cramersv, as.factor(vt.cv0$sRc), p.adjust = "bonf"))
