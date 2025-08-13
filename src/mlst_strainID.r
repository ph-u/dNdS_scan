#!/bin/env Rscript
# author: ph-u
# script: mlst_strainID.r
# desc: Summarize which isolates are likely the same strain
# in: Rscript mlst_strainID.r
# out: ../res/mlst_strainID.{pdf,csv}
# arg: 0
# date: 20250412

source("p_src.r")
set.seed(123)
mlst.st = read.csv(paste0(pT[2],"mlst_realphy--metadata.csv"), header = T)
f0 = list.files(pT[1],"mlst_stSNP", full.names = T)
st.uniq = unique(mlst.st$seqType)

cat(date(),": Get gene lengths (mlst_strainID.r)\n")
gLen = rep(NA,nrow(f));for(i in 1:nrow(f)){
  cat(date(),": Reading",i,"/",nrow(f),";",round(i/nrow(f)*100,2),"%     \r")
  gLen[i] = length(as.character(read.FASTA(paste0(pT[3],f$fNam[i]),type = "DNA"))[[1]])}
write.csv(data.frame(gNam = f$gEne, gLen = gLen), paste0(pT[1],"mlst_strainID--gLen.csv"), row.names = F, quote = F)

cat(date(),": Start SNP plotting (mlst_strainID.r)\n")
pdf(paste0(pT[2],"mlst_strainID.pdf"), width = 10, height = 20, paper = "a4")
par(mar = c(5,4,4,0)+.1, mfrow = c(4,1))
i5 = 1;for(i in 1:length(st.uniq)){
  cat(date(),": Standardizing",paste0("ST",st.uniq[i]),";",i,"/",length(st.uniq),"\n")
  i0 = read.csv(f0[grep(paste0("ST",st.uniq[i]),f0)], header = T, row.names = 1)/gLen
  i0$gNam = GTF$gNam[match(row.names(i0),GTF$Locus.Tag)]
  i0$gNam[is.na(i0$gNam)] = row.names(i0)[is.na(i0$gNam)]
  cat(date(),": Plotting",paste0("ST",st.uniq[i]),";",i,"/",length(st.uniq),"\n")
#  i1 = colSums(i0, na.rm = T)
#  matplot(i0*100, type = "l", ylab = "SNP percentage", xlab = "Gene Position", pch = 1:ncol(i0), col = paste0(substr(cBp,1,7),"44"), main = paste0("ST",st.uniq[i]), lty = 1)
  i1 = mlst.st[which(mlst.st$seqType==st.uniq[i]),1]
  i4 = c();for(i2 in 1:(length(i1)-1)){for(i3 in (i2+1):length(i1)){i4 = c(i4,paste0(i1[i2]," vs ",i1[i3],";\n(",paste0(mlst.st[which(mlst.st[,1]==i1[i2]),c(2,5)], collapse = ","),") vs (",paste0(mlst.st[which(mlst.st[,1]==i1[i3]),c(2,5)], collapse = ","),")"))}}
  for(i1 in 1:length(i4)){
    plot(x = 1:nrow(i0), y = i0[,i1]*100, main = paste0("ST",st.uniq[i],": ",i4[i1],": col ",i1,": G",i5), type = "p", pch = 19, ylab = "SNP percentage", xlab = "Gene Position", ylim = c(0,100))
    i6 = ifelse(i0[,i1]<.8,"",i0$gNam)
    if(i5==1){
      r0 = i6
    }else if(i5==2){
      r0 = data.frame(G1=r0,G2=i6)
      row.names(r0) = row.names(i0)
    }else{
      r0[,i5] = i6
    };i5 = i5 + 1
  }
#  legend("topright", legend = i4, fill = cBp, ncol = 3)
};rm(i,i0,i1,i2,i3,i4,i5,i6)
invisible(dev.off())
colnames(r0) = paste0("G",1:ncol(r0))
write.csv(r0[which(apply(r0,1,function(x){!all(x=="")})==T),], paste0(pT[2],"mlst_strainID.csv"), quote = F, row.names = T)
cat(date(),": mlst_strainID.r Done\n")
