#!/bin/env Rscript
# author: ph-u
# script: dNdS_reCon_srcDiff-str.r
# desc: summarize reCon result with sampling source comparison
# in: Rscript dNdS_reCon_srcDiff-str.r [PANum] [statistics_colname]
# out: data/dNdS_[PANum]_reCon_srcDiff-str.csv
# arg: 2
# date: 20240724

argv = (commandArgs(T))
#argv = c("PA2934","mean")
source("p_src.r")

rDNDS = read.csv(paste0(pT[3],list.files(pT[3],"-rDNDS")[grep(paste0("_",argv[1],"-"), list.files(pT[3],"-rDNDS"))]), header = T)
rDNDS$dNdS[is.infinite(rDNDS$dNdS)] = max(rDNDS$dNdS[is.finite(rDNDS$dNdS)])+2
reCon = read.csv(paste0(pT[3],list.files(pT[3],"-reCon")[grep(paste0("_",argv[1],"-"), list.files(pT[3],"-reCon"))]), header = T)

##### Reconstruction of strain level residue-resolved dNdS #####
cat("Start reconstructing residue-level dNdS on each strain:",date(),"\n")
for(i in 1:nrow(mEta)){
    cat(i,"/",nrow(mEta),":",date(),"     \r")
    rD = rDNDS[grep(mEta$assemblyInfo.genbankAssmAccession[i], rDNDS$clinical),]
    r.t = reCon[,-(1:3)]
    colnames(r.t) = gsub("dNdS", mEta$assemblyInfo.genbankAssmAccession[i], colnames(r.t))
    for(i0 in 1:nrow(r.t)){ r.t[i0,] = summary(rD$dNdS[which(rD$ntStart<=reCon$ntPos[i0] & rD$ntEnd>=reCon$ntPos[i0] & !is.na(rD$dNdS))]) }
    if(i>1){r0 = cbind(r0, r.t)}else{r0 = cbind(reCon[,1:3],r.t)}
};rm(i, i0, rD, r.t)
cat("Done reconstruction:",date(),"\n")

write.csv(r0, paste0(pT[1],"dNdS_reCon_",argv[1],"_str.csv"), row.names = T, quote = F)

##### statistical comparison between CF & env #####
d.cat.c = unique(mEta$sOurce)[order(unique(mEta$sOurce))]
d.cat = list(mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==d.cat.c[1])], mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==d.cat.c[2])])
names(d.cat) = d.cat.c[1:2]

r1.c = c("diff", "wilcox", "p.val")
r1 = as.data.frame(matrix(nr = nrow(reCon), nc = length(r1.c)))
colnames(r1) = r1.c
r1.raw = r0[,grep(argv[2], colnames(r0))]
colnames(r1.raw) = gsub(paste0(".",argv[2]),"",colnames(r1.raw))
for(i in 1:nrow(r1)){
    r1.r = list(as.numeric(r1.raw[i,which(colnames(r1.raw) %in% d.cat[[1]])]), as.numeric(r1.raw[i,which(colnames(r1.raw) %in% d.cat[[2]])]))
    r1.t = summary(r1.r[[2]])-summary(r1.r[[1]])
    r1.s = wilcox.test(r1.r[[1]], r1.r[[2]])
    r1[i,] = c(r1.t[which(tolower(names(r1.t))==argv[2])], r1.s$statistic, r1.s$p.value)
};rm(i, r1.r, r1.s, r1.t)
r1$p.adj = p.adjust(r1$p.val, method = "BH")
r1$absLOG10 = -log10(r1$p.adj)*ifelse(r1$diff<0, -1,1)
#r1$sCaled = (r1$absLOG10-min(r1$absLOG10))/max(r1$absLOG10)

write.csv(cbind(reCon[,1:3],r1), paste0(pT[1],"dNdS_reCon_",argv[1],"_strP.csv"), row.names = F, quote = F)

##### Manhattan plot #####
aCtive = c(129, 153, 297)
pdf(paste0(pT[2],"dNdS_reCon_",argv[1],"_strP.pdf"), width = 5, height = 5)
par(mar = c(5,4,0,1)+.1)
plot(x = 1:nrow(r1), y = r1$absLOG10, type = "l", lwd = 3, xlab = paste0(argv[1]," residue position"), ylab = bquote(~-Log[10] ~ italic(P.adj)), yaxt = "n", col = cBp[1], ylim = c(-2.5, 2))
points(x = aCtive, y = r1$absLOG10[aCtive])
axis(2, at = -2:2, labels = c(2,1,0,1,2))
abline(h = c(-1,1), lty = 3)
rect(xleft = c(25,155), ybottom = -2.5, xright = c(319,242), ytop = -2.3, lty = c(2,1))
text(x = c(80, 200, 270), y = -2.4, labels = c("core", "cap", "core"))
text(x = c(129, 153, 297)+10, y = r1$absLOG10[aCtive]+.2, labels = c("D129", "E153", "H297"))
mtext(d.cat.c[1:2], 4, at = c(1.5,-1.5))
invisible(dev.off())
