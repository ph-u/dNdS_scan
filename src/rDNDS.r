#!/bin/env Rscript
# author: ph-u
# script: rDNDS.r
# desc: calculate & plot residue-wise dN/dS contribution of the given gene
# in: Rscript rDNDS.r [srcFile] [prefix] [gene] [locus_Tag]
# out: res/02_[prefix]_[gene]_[locus_Tag].{csv,pdf}
# arg: 4
# date: 20240119

argv=(commandArgs(T))
#argv = c("../data/01_PAO1_107_PA0506.csv","PAO1_107","PA0506","fadE1")
aRg = gsub("01_","00_",argv[1])
library(ape); library(msa)
rAin = rev(rainbow(130)[1:100])

a0 = read.csv(argv[1], header = T)
a0 = a0[which(!is.na(a0$dNdS)),]

##### Modify substitution matrix for conservation score #####
data(BLOSUM62) # https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/Biostrings/html/substitution_matrices.html
BLOSUMxx = matrix(numeric(prod(dim(BLOSUM62))), nr = nrow(BLOSUM62), nc = ncol(BLOSUM62))
row.names(BLOSUMxx) = colnames(BLOSUMxx) = row.names(BLOSUM62)
for(i in 1:nrow(BLOSUMxx)){BLOSUMxx[i,i] = 1};rm(i)

##### Segregate CF & non-CF samples #####
#mD = read.csv("../../1_01_PA2668/data/metaData.csv", header=T)
#mD[mD=="missing"]=mD[mD==""]=NA

##### f: get rows for necessary category change #####
#cHg = function(tErms, tArget, df = mD){ 
#        i0 = c() 
#        for(i in tErms){i0 = unique(c(i0,grep(i,df[,grep("isolation_", colnames(df))])))};rm(i)
#        df$sOurce[i0] = tArget
#        return(df)}

#mD$cOuntry = mD$sOurce = NA
#for(i in 1:nrow(mD)){mD$cOuntry[i] = ifelse(is.na(mD[i,grep("geo_loc", colnames(mD))]),"Unknown",strsplit(mD[i,grep("geo_loc", colnames(mD))],":")[[1]][1])};rm(i)
#mD = cHg(c("ystic ", "sputum"), "Cystic fibrosis")
#mD = cHg(c("skin", "toe"), "Skin")
#mD = cHg(c("iver", "water", "soil", "milk", "nvironm", "egg"), "Environmental")
#mD$sOurce[which(is.na(mD[,grep("isolation_", colnames(mD))]))] = "Unknown"
#mD$sOurce[is.na(mD$sOurce)] = "Other infections"

##### Get AA seq alignment #####
if(!file.exists(gsub(".csv$", "-aaAlign.fa", aRg))){
    writeXStringSet(unmasked(msa(readAAStringSet(gsub(".csv$","_dbAA.fa", aRg)), "ClustalOmega")), file=gsub(".csv$", "-aaAlign.fa", aRg))
}

##### Reconstruction residue-wise dN/dS & protein conservation #####
iX = "All";tG = unique(a0$clinical); a = a0
#i2 = c("Cystic fibrosis", NA)
#for(i1 in 1:length(i2)){
#    if(!is.na(i2[i1])){
#        iX = gsub(" ","",i2[i1])
#        tG = mD$assemblyInfo.genbankAssmAccession[which(mD$sOurce==i2[i1])]
#    }else{
#        iX = "Others"
#        tG = c(); for(i3 in i2[-i1]){tG = c(tG, mD$assemblyInfo.genbankAssmAccession[which(mD$sOurce==i3)])};rm(i3)
#    }
#    cat("Extracting samples:",iX,date(),"\n")
#    xMD = c(); for(i3 in 1:length(tG)){xMD = c(xMD, grep(tG[i3], a0$clinical))};rm(i3)
#    a = a0[xMD,]

##### Get AA alignment of selected sample source #####
    cat("Extracting peptides from source:",iX,date(),"\n")
    aL0 = as.character(read.FASTA(gsub(".csv$", "-aaAlign.fa", aRg), type="AA"))
    aLn = vector(mode="list", length = length(tG))
    for(i3 in 1:length(tG)){
        x1 = grep(tG[i3], names(aL0))
        aLn[[i3]] = aL0[[x1]]
        names(aLn)[i3] = names(aL0)[x1]
    };rm(i3, x1)
    write.FASTA(as.AAbin(aLn), gsub(".csv$",paste0("-aaAlign-",iX,".fa"),aRg))
    aLn = msa(readAAStringSet(gsub(".csv$",paste0("-aaAlign-",iX,".fa"), aRg)), "ClustalOmega")
    writeXStringSet(unmasked(aLn), file=gsub(".csv$",paste0("-aaAlign-",iX,".fa"),aRg))
    aL0 = msaConservationScore(aLn, BLOSUMxx)

# ranged from 0 (all dS) to inf (all dN)
    cat("Summarizing dN/dS from designated samples:",iX,date(),"\n")
    aMino = data.frame(ntPos=seq(1, max(c(a$start,a$end)), 3), mEdian=NA, nRef=NA, p.val=NA)

##### max number of available dNdS values #####
    sqMax = 1:nrow(aMino)
    sqMax[(ceiling(nrow(aMino)/2)+!(nrow(aMino)%%2)):nrow(aMino)] = ceiling(nrow(aMino)/2):1
    sqMax[sqMax>(a$end[1]-a$start[1])] = (a$end[1]-a$start[1])
    sqMax = sqMax*854 #length(list.files(".", ".ndb"))

##### Protein residue essentiality potential #####
    cat("Reconstructing residue-wise dN/dS:",iX,date(),"\n")
    for(i in 1:nrow(aMino)){
        i0 = which(a$start <= aMino$ntPos[i] & a$end >= aMino$ntPos[i])
        aMino$mEdian[i] = median(a$dNdS[i0])
        aMino$nRef[i] = length(i0)
        aMino$p.val = wilcox.test(a$dNdS[i0], a$dNdS)$p.value
    };rm(i, i0)
    aMino$p.adj = p.adjust(aMino$p.val, method = "BH")
    aMino$mEdian[which(aMino$p.adj > .05)] = median(a$dNdS)
    aMino$cOl = gray(1-(.1+.9*aMino$nRef/sqMax))
    aMino$prot = c(1-aL0/max(aL0),0) # last codon is stop codon, no animo acid but conserved

    pdf(paste0("../res/02_",argv[2],"_",argv[3],"_",argv[4],"_",iX,".pdf"), width = 150)
    plot(0,0,col="#00000000", xlab = "Residue position", ylab = "Median dN/dS Contribution", main = paste(argv[3],"/",argv[4]), xlim = c(1,nrow(aMino)), ylim = c(0, max(1,aMino$mEdian)+.2))
    for(i in 1:100){
        rect(0, diff(range(aMino$mEdian))*(i-1)/100, nrow(aMino), diff(range(aMino$mEdian))*i/100, col = rAin[i], border = NA)
    };rm(i)
    points(1:nrow(aMino), aMino$mEdian, col=aMino$cOl, type="p", pch=19)
    lines(x=1:nrow(aMino), y=aMino$mEdian, col="#000000ff", lty=1)
    text(x = 1:nrow(aMino), y = max(aMino$mEdian)+.15, cex = 2, labels = strsplit(msaConsensusSequence(aLn),"")[[1]])
    abline(h=1, col="#ffffffff", lty=2, lwd=3)
    invisible(dev.off())

    write.csv(aMino, paste0("../res/02_",argv[2],"_",argv[3],"_",argv[4],"_",iX,".csv"), row.names = F, quote = F)
#}; rm(i1)
