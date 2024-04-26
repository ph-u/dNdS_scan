#!/bin/env Rscript
# author: ph-u
# script: varTypeBar.r
# desc: sequence variation type in all PAO1 protein-coding genes
# in: Rscript varTypeBar.r
# out: res/vt_Bar.pdf
# arg: 0
# date: 20240320

source("p_src.r")
f = f[which(f$gEne %in% protList),]

cat("Categorizing variation types:",date(),"\n")
for(i0 in 1:nrow(f)){
    cat("Processing",f$gEne[i0],":",date(),"     \r")
    eSx = unique(read.csv(paste0(pT[3],paste0("PAO1_107_",f$gEne[i0],"--dbSum.csv")), header = T))
    eSx$gene = f$gEne[i0]
    vType = strsplit(eSx$varType, ":")
    eSx$vType = eSx$varType
    for(i in 1:length(vType)){if(length(vType[[i]]) > 1){eSx$vType[i] = paste(vType[[i]][c(1,length(vType[[i]]))], collapse = ":")}};rm(i, vType)
    if(i0==1){eS = eSx}else{eS = rbind(eS,eSx)}
};rm(i0,eSx)

cat("\nPlotting:",date(),"\n")
x = length(unique(eS$clinical)) # ceiling(x/(10^nchar(x)))*10^nchar(x)
x0 = 8; x1 = round(nrow(f)/x0)
pdf(paste0(pT[2],"vt_Bar.pdf"), width = x1, height = x0*7)
par(mar = c(5,4,1,10)+.1, mfrow = c(x0,1), xpd = T)
teS = table(eS[,c("vType","gene")])/x*100
for(i in 1:x0){
    barplot(teS[,(x1*(i-1)+1):min(x1*i,nrow(f))], col = cBp, ylab = "IPCD isolates (%; 100% = 854 isolates)")
};rm(i)
legend("topleft", legend = row.names(teS), pch = 20, cex = 6, col = cBp, inset = c(0,0))
invisible(dev.off())
cat("Exported barplot:",date(),"\n")

##### List top 20 genes of each category #####
#for(i in 1:nrow(teS)){
#    cat("Top 20 genes",row.names(teS)[i],"to ref seq: ", paste(colnames(teS)[rev(order(teS[i,]))[1:20]], collapse = ", "), "\n")
#};rm(i)
write.csv(t(teS), paste0(pT[2],"vt_Bar.csv"), row.names = T, quote = F)
