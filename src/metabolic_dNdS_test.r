#!/bin/env Rscript
# author: ph-u
# script: metabolic_dNdS_test.r
# desc: Does the living environment of PA affect the evolutionary potential of genes contributing to one metabolic pathway?
# in: Rscript metabolic_dNdS_test.r
# out: NA
# arg: 0
# date: 20240403

p0 = 1

source("p_metabolism_PAO1.r")
library(PMCMRplus)
gBioc = gBioc[!(gBioc$metabolism %in% c("Others","Multiple")),]

if(p0!=1){
cat("Mapping dN/dS values:", date(), "\n")
for(i in 1:nrow(gBioc)){
    cat(i, "/", nrow(gBioc), "(", round(i/nrow(gBioc)*100), "% ;", gBioc$protList[i], ")", date(), "     \r")
    i0 = read.csv(paste0(pT[3], "PAO1_107_", gBioc$protList[i], "--dbSum.csv"), header = T)
    if(i==1){
        r0 = as.data.frame(matrix(nr = nrow(i0), nc = nrow(gBioc)+1))
        row.names(r0) = read.table(text = gsub("_ASM", "@", i0$clinical), sep = "@")[,1]
        colnames(r0) = c(gBioc$protList,"sRc")
    }
    row.names(i0) = read.table(text = gsub("_ASM", "@", i0$clinical), sep = "@")[,1]
    r0[,gBioc$protList[i]] = i0$dNdS[match(row.names(i0), row.names(r0))]
};rm(i); cat("\nMapping done:", date(), "\n")
r0$sRc = mEta$sOurce[match(mEta$assemblyInfo.genbankAssmAccession, row.names(r0))]

cat("Rearrange dN/dS values:", date(), "\n")
for(i in 1:(ncol(r0)-1)){
    cat(i, "/", ncol(r0)-1, "(", round(i/(ncol(r0)-1)*100), "% ;", colnames(r0)[i], ")", date(), "     \r")
    r1.0 = r0[,c(i,ncol(r0))]
    r1.0$gene = colnames(r1.0)[1]
    colnames(r1.0)[1] = "dNdS"
    r1.0$clinical = row.names(r1.0)
    r1.0$metabolism = gBioc$metabolism[i]
    if(i==1){ r1 = r1.0 }else{ r1 = rbind(r1, r1.0) }
};rm(i, r1.0); cat("\nRearrangemnet done:", date(), "\n")
}else{
    r0 = read.csv(paste0(pT[1], "metabolic_dNdS_test_wide.csv"), header = T, row.names = 1)
    r1 = read.csv(paste0(pT[1], "metabolic_dNdS_test.csv"), header = T)}

##### Statistics test: resolved #####
pairwise.wilcox.test(r1$dNdS, r1$sRc, p.adjust = "BH")
## do not know whether samples are independent, so KW is inappropriate (https://stats.stackexchange.com/questions/302464/wilcoxon-signed-rank-vs-kruskal-wallis)
#print(kruskal.test(r1$dNdS ~ r1$sRc))
#print(kW <- kwAllPairsNemenyiTest(r1$dNdS ~ as.factor(r1$sRc)))
#kWp = data.frame(g1 = rep(colnames(kW$p.value), each = nrow(kW$p.value)), g2 = rep(row.names(kW$p.value), nrow(kW$p.value)), p.adj = p.adjust(kW$p.value, method = "BH"))
if(p0!=1){
    write.csv(r0, paste0(pT[1], "metabolic_dNdS_test_wide.csv"), row.names = T, quote = F)
    write.csv(r1, paste0(pT[1], "metabolic_dNdS_test.csv"), row.names = F, quote = F)}

##### Distribution plot & summaries #####
x0 = unique(r1$sRc)
x0.1 = paste0("dNdS.", c("min", "1Q", "2Q", "mean", "3Q", "max", "NA"))
x1 = as.data.frame(matrix(nr = length(x0), nc = length(x0.1)))
row.names(x1) = x0; colnames(x1) = x0.1
pdf(paste0(pT[2], "metabolic_dNdS_test.pdf"))
for(i in 1:length(x0)){
    x0.1 = r1$dNdS[grep(x0[i], r1$sRc)]
    x1[i,] = summary(x0.1, na.rm = T)
    if(i==1){
        hist(x0.1, col = sub("FF$","77",cBp[i]), freq = T, xlab = "dN/dS")
    }else{
        hist(x0.1, col = sub("FF$","77",cBp[i]), freq = T, add = T)
}};rm(i, x0.1)
legend("topright", legend = x0, pch = 19, col = sub("FF$","77",cBp[1:length(x0)]))
invisible(dev.off())

pdf(paste0(pT[2], "metabolic_dNdS_test_box.pdf"), width = length(x0)*1.6, height = 10)
par(mar = c(4,4,0,0)+.1, mfrow = c(2,1))
boxplot(r1$dNdS ~ as.factor(r1$sRc), col = "#00000000", xlab = "", ylab = "dN/dS")
boxplot(r1$dNdS ~ as.factor(r1$sRc), col = "#00000000", xlab = "", ylab = "dN/dS", ylim = c(0,.15))
invisible(dev.off())

##### PCA with NA on dN/dS by sample source #####
#r0.b = r0
#x2 = unique(r1$metabolism)
#for(i in 1:length(x2)){
#    cat(i,"/",length(x2),"(",round(i/length(x2)*100),"% ;",x2[i],")",date(),"                   \r")
#    x1 = r0[,c(colnames(r0)[which(colnames(r0) %in% unique(r1$gene[which(r1$metabolism==x2[i])]))],"sRc")]
## clean inf values (https://stackoverflow.com/questions/12188509/cleaning-inf-values-from-an-r-dataframe)
#    x1 = do.call(data.frame,lapply(x1, function(x) replace(x, is.infinite(x),max(x1[,-ncol(x1)], na.rm = T)+1)))
#    rRM = c();for(i0 in 1:(ncol(x1)-1)){if(length(unique(x1[,i0]))<3){rRM = c(rRM,i0)}};rm(i0)
#    if(length(rRM)>0){x1 = x1[,-rRM]}
#    rRM = c();for(i0 in 1:nrow(x1)){if(sum(is.na(x1[i0,]))==(ncol(x1)-1)){rRM = c(rRM,i0)}};rm(i0)
#    if(length(rRM)>0){x1 = x1[-rRM,]}

#    if(ncol(x1)>2){
#        pPca = pca(as.matrix(x1[,-ncol(x1)]), method = "ppca")
#        pPlt = pcaSwap(as.matrix(x1[,-ncol(x1)]))
#        p0 = ggbiplot(pPlt, obs.scale = 1, var.scale = 1, groups = as.factor(x1$sRc), ellipse = TRUE, ellipse.prob = 0.95) + theme_bw() +
#            scale_color_manual(values=setNames(cBp[1:length(x0)], x0)) +
#            guides(color=guide_legend(title="Sample source")) +
#            theme(legend.position = 'bottom', legend.direction = "vertical",
#                panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                legend.title = element_text(size=18),
#                axis.text=element_text(size=16),axis.title=element_text(size=14),
#                plot.margin=margin(t=0,r=0,b=0,l=1))
#        ggsave(paste0(pT[2],"metabolism_pca/mb_pca_",gsub(" ","-",x2[i]),".pdf"), plot = pcaLAB(p0,round(pPca@R2*100,1)), width = 6, height = 7)
#}else{cat("Too few member genes, metabolism skipped:",x2[i],"-",date(),"\n")}};rm(i)
