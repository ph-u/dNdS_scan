#!/bin/env Rscript
# author: ph-u
# script: goMap_wCompare.r
# desc: compare whole gene dN/dS for ones with genomegaMap result
# in: Rscript goMap_wCompare.r
# out: res/gmCompare/gm_*--box.pdf
# arg: 0
# date: 20240515

source("p_src.r")
f.goMap = list.files(pT[3],"_whole.csv")
f.rm = c();for(i in 1:length(f.goMap)){if(file.size(paste0(pT[3],f.goMap[i]))==0){f.rm = c(f.rm, i)}};f.goMap = f.goMap[-f.rm];rm(i, f.rm)
f.dbSum = list.files(pT[3],"--dbSum")
f.paNum = read.table(text = f.goMap, sep = "_")[,4]

##### Collect dN/dS from genomegaMap & program #####
r0.c = c("PAnum", "W", "p.val")
r0 = as.data.frame(matrix(nr = length(f.goMap), nc = length(r0.c)))
colnames(r0) = r0.c
r0[,1] = f.paNum

cat("Analyzing dN/dS differences:", date(), "\n")
for(i in 1:length(f.goMap)){
    cat("Scanning",f.paNum[i], "(", i, "/", length(f.goMap), ";", round(i/length(f.goMap)*100), "% ) :", date(), "     \r")
    f.gm = read.table(paste0(pT[3], f.goMap[i]), header = T, sep = "\t")
    f.db = read.csv(paste0(pT[3], f.dbSum[grep(paste0(f.paNum[i],"-"), f.dbSum)]), header = T)
    f.gm = f.gm[which(f.gm[,2] > max(f.gm[,2])-1),]#[rev(order(f.gm[,2])),][1:nrow(f.db),] #[which(f.gm[,2] > max(f.gm[,2])-.7),]

##### Assumptions #####
    f.db$dNdS[which(f.db$varType=="identical")] = 0
#    f.db = f.db[!is.na(f.db$dNdS),]
    f.gd = data.frame(src=c(rep("genomegaMap", nrow(f.gm)), rep("IPCD", nrow(f.db))), dNdS = c(f.gm$omega0,f.db$dNdS))
    if(all(is.na(f.db$dNdS))){kW = list(statistic=NA, p.value=NA)}else{
        kW = wilcox.test(f.gd$dNdS[which(f.gd$src=="IPCD")], f.gd$dNdS[which(f.gd$src!="IPCD")])
    }
    f.gd$dNdS[is.infinite(f.gd$dNdS)]=max(f.gd$dNdS, na.rm=T, 1)+.1
    r0[i,-1] = c(unname(kW$statistic,0), kW$p.value)

    pdf(paste0(pT[2],"gmCompare/gm_",f.paNum[i],"--box.pdf"))
    boxplot(f.gd$dNdS~f.gd$src, pch = 4, col = "#00000000", xlab = paste("Data Source :",f.paNum[i]), ylab = "dN/dS")
    segments(x0 = 1, y0 = max(f.gd$dNdS, na.rm=T)*.8, x1 = 2, lwd = 3)
    text(1.5, max(f.gd$dNdS, na.rm=T)*.82, label=paste("Wilcox W =", unname(round(kW$statistic,0)), "; p.value", ifelse(kW$p.value < .001, "<< 0.01", ifelse(kW$p.value > .05, "> 0.05", ifelse(kW$p.value < .01, "< 0.01", round(kW$p.value,2))))))
    invisible(dev.off())
};rm(i); cat("\nDone scanning:", date(), "\n")
write.csv(r0, paste0(pT[2],"gm_Compare.csv"), row.names = F, quote = F)
