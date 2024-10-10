#!/bin/env Rscript
# author: ph-u
# script: clinical_gSeq_Map.r
# desc: map linearly gene neighbourhood of IPCD initial collection
# in: Rscript clinical_gSeq_Map.r
# out: res/clinical_gSeq_Map.pdf
# arg: 0
# date: 20240823

rUn = 1
source("p_src.r")
dbSum = list.files(pT[3],"--dbSum")

##### Collect seq start positions #####
rSeq = seq(1,2*length(dbSum),2)
if(rUn == 0){
    cat("Collecting locus and respective sequence start positions:",date(),"\n")
    for(i in 1:length(dbSum)){
        cat("Processing gene",f$gEne[i],":",date(),"     \r")
        i.0 = as.data.frame(t(read.csv(paste0(pT[3],dbSum[i]), header = T)[,c(1,2,4)]))
        colnames(i.0) = i.0[1,]
        i.0 = i.0[,order(colnames(i.0))]
        if(i==1){r0 = as.data.frame(matrix(nr = 2*length(dbSum), nc = ncol(i.0))); colnames(r0) = colnames(i.0)}
        r0[rSeq[i]+(0:1),] = i.0[-1,]
    };rm(i,i.0);cat("\nCollection completed:",date(),"\n")
    row.names(r0) = paste0(rep(f$gEne, each = 2),".",c("locus","start"))
    write.csv(r0,paste0(pT[1],"clinical_gSeq_Map.csv"), quote = F, row.names = T)
}else{
    cat("Reading locus and respective sequence start positions:",date(),"\n")
    r0 = read.csv(paste0(pT[1],"clinical_gSeq_Map.csv"),header = T, row.names = 1)
}

##### Strain genome mapping #####
gMap = as.data.frame(matrix(nr = length(rSeq), nc = ncol(r0)))
colnames(gMap) = colnames(r0)
row.names(gMap) = f$gEne
cat("Collecting seq positions map coordinates:",date(),"\n")
for(i in 1:ncol(gMap)){
    cat("Processing strain",colnames(gMap)[i],"(",round(i/ncol(gMap)*100,1),"% ):",date(),"     \r")
    gLocus = as.numeric(as.factor(r0[rSeq,i]))
    gLocus[is.na(gLocus)] = 0
    gStart = order(r0[rSeq+1,i])
    gMap[,i] = gLocus+gStart/(10^ceiling(log10(length(gStart))))
};rm(i);cat("\nCollection completed:",date(),"\n")
write.csv(gMap,paste0(pT[1],"clinical_gMap.csv"), quote = F, row.names = T)

##### Genome mapping #####
colnames(gMap) = read.table(text = gsub("_ASM","@",colnames(gMap)), sep = "@")[,1]

#i.0 = mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==unique(mEta$sOurce)[grep("Cystic",unique(mEta$sOurce))])]
i.0 = c();for(i in unique(mEta$sOurce)[order(unique(mEta$sOurce))]){i.0 = c(i.0,mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==i)])};rm(i)

#gMap = gMap[,which(colnames(gMap) %in% i.0)]
gMap = gMap[,order(match(colnames(gMap), i.0))]

cat("Graphics visual-processing:",date(),"\n")
for(i in 1:ncol(gMap)){
    cat("Processing strain",colnames(gMap)[i],"(",round(i/ncol(gMap)*100,1),"% ):",date(),"     \r")
    gMap[,i] = ifelse(gMap[,i]<1,gMap[,i],gMap[,i]/max(gMap[,i])*ceiling(max(gMap)))
};rm(i);cat("\nProcessing completed:",date(),"\n")

PAO1.gPos = (1:nrow(gMap))/nrow(gMap)*ceiling(max(gMap))
cHigh = grep("PA2125|PA2384",f$gEne)

for(i0 in 1:2){
    if(i0==1){
        cOl.gPos = paste0(rainbow(nrow(gMap)),"33", sep = "")
        cOl.gPos[cHigh[1]:cHigh[2]] = "#000000ff"
    }else{
        cOl.gPos = rep("#00000033", nrow(gMap))
        cOl.gPos[cHigh[1]:cHigh[2]] = paste0(rainbow(diff(cHigh)+1),"ff", sep = "")
    }

    #pdf(paste0(pT[2],"clinical_gSeq_Map.pdf"), height = floor(ncol(gMap)/8)*2, width = 15)
    pdf(paste0(pT[2],"clinical_gSeq_Map--all",i0,".pdf"), height = floor(ncol(gMap)/8)*2, width = 15)
    par(mar = c(5,15,1,0)+.1)
    plot(c(0,ceiling(max(gMap))),c(0,ncol(gMap)+1), col = "#00000000", xlab = "Scaffold-gene order", ylab = "", yaxt = "n")
    axis(2, at = 0:ncol(gMap), labels = c("PAO1-ref",paste0(colnames(gMap)," - ",mEta$sOurce[match(colnames(gMap),mEta$assemblyInfo.genbankAssmAccession)])), las=1, cex = .5)
    abline(v = 1, col = cBp[1])

    cat("Plotting seq positions map coordinates (scaffold background):",date(),"\n")
    for(i in 1:nrow(gMap)){ # 1:nrow(gMap) 2100:2450
        cat("Processing gene",f$gEne[i],"(",round(i/nrow(gMap)*100,1),"% ):",date(),"     \r")
        points(x = c(PAO1.gPos[i],gMap[i,]), y = 0:ncol(gMap), pch = 3, cex = .1, col = cOl.gPos[i])
    };rm(i);cat("\nPlot background done:",date(),"\n")

    cat("Highlight region of interest:",date(),"\n")
    for(i in cHigh[1]:cHigh[2]){ # 1:nrow(gMap) 2100:2450
        cat("Processing gene",f$gEne[i],"(",round((i-cHigh[1]+1)/(diff(cHigh)+1)*100,1),"% ):",date(),"     \r")
        points(x = c(PAO1.gPos[i],gMap[i,]), y = 0:ncol(gMap), pch = 3, cex = .1, col = cOl.gPos[i])
    };rm(i);cat("\nPlot ready to be exported:",date(),"\n")

    invisible(dev.off())
};rm(i0)
cat("Plot exported:",date(),"\n")
