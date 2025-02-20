#!/bin/env Rscript
# author: ph-u
# script: dNdS_infect-vol.r
# desc: volcano plot of infectious vs non-infectious dN/dS
# in: Rscript dNdS_infect-vol.r
# out: data/dNdS_infect-vol.csv, res/dNdS_volcano_[cond1]_[cond2]-vol.pdf
# arg: 0
# date: 20240706

source("p_src.r")
d.nam = paste0(pT[1],"dNdS_infect-vol.csv")
d.c = unique(mEta$sOurce)

library(EnhancedVolcano)

if(!file.exists(d.nam) | file.info(d.nam)$size < 3124638){
##### f: dN/dS volcano plot equation #####
    f.dNdS.volcano = function(x){
        x0 = c(log10(median(x[[1]]+1e-16, na.rm=T)), log10(median(x[[2]]+1e-16, na.rm=T)))
        return(ifelse(x0[1]<x0[2],-1,1)*abs(log2(abs(x0[1]/x0[2]))))} # condition 1 more important (tend towards dN/dS conserved values) on the left, condition 2 more important on the rihgt

##### Collect log2[dNdS-median]C & statistics from all genes #####
    f0 = list.files(pT[3], "dbSum")
    f = data.frame(fNam = f0, gEne = read.table(text = gsub("_PA","@PA",gsub("--","@",f0)), sep = "@")[,2])

    d.cat = vector(mode = "list", length = length(d.c))
    names(d.cat) = d.c
    for(i in 1:length(d.cat)){d.cat[[i]] = mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==names(d.cat)[i])]};rm(i)
    d.Comp = list(c(d.cat[[2]],d.cat[[4]]), c(d.cat[[3]],d.cat[[5]])) # infectious vs non

    cat("Starting dNdS statistics data-collection:",date(),"\n")
    r0.c = c("gene", "cond1", "cond2", "log2FC", "p.val")
    r0 = as.data.frame(matrix(nr = choose(length(d.c)-1,2)+1, nc = length(r0.c)))
    colnames(r0) = r0.c

    for(i in 1:nrow(f)){
        cat("Processing gene",f$gEne[i],"(",i,"/",nrow(f),") :",date(),"       \r")
        d0 = read.csv(paste0(pT[3], f$fNam[i]), header = T)
        d0$dNdS[is.infinite(d0$dNdS)] = Inf
        d0$clinical = read.table(text = gsub("_ASM","@",d0$clinical), sep = "@")[,1]

        r0[,1] = f$gEne[i]
        r0[1,2:3] = c("infection","non-infection")
        r0.L = list(d0$dNdS[which(d0$clinical %in% d.Comp[[1]])], d0$dNdS[which(d0$clinical %in% d.Comp[[2]])])
        if(sum(is.na(r0.L[[1]]))<length(r0.L[[1]]) & sum(is.na(r0.L[[2]]))<length(r0.L[[2]])){
            r0[1,4:5] = c(f.dNdS.volcano(r0.L), wilcox.test(r0.L[[1]], r0.L[[2]])$p.value)
        }else{r0[1,4:5] = rep(NA,2)}
        i2 = 2;for(i0 in 2:(length(d.cat)-1)){ for(i1 in (i0+1):length(d.cat)){
            r0[i2,2:3] = names(d.cat)[c(i0,i1)]
            r0.L = list(d0$dNdS[which(d0$clinical %in% d.cat[[i0]])], d0$dNdS[which(d0$clinical %in% d.cat[[i1]])])
            if(sum(is.na(r0.L[[1]]))<length(r0.L[[1]]) & sum(is.na(r0.L[[2]]))<length(r0.L[[2]])){
                r0[i2,4:5] = c(f.dNdS.volcano(r0.L), wilcox.test(r0.L[[1]], r0.L[[2]])$p.value)
            }else{r0[1,4:5] = rep(NA,2)}
            i2 = i2+1
        }}
        if(i==1){write.csv(r0, d.nam, quote = F, row.names = F)}else{write.csv(rbind(read.csv(d.nam, header = T),r0), d.nam, quote = F, row.names = F)}
    };rm(i, i0, i1, d0, r0.L)

## p.adj calculation
    cat("\np.adj calculation:",date(),"\n")
    r0 = read.csv(d.nam, header = T)
    r0$cOnd = paste0(r0$cond1,r0$cond2)
    r0.c = unique(r0$cOnd)
    r0$p.adj = NA
    for(i in 1:length(r0.c)){r0$p.adj[which(r0$cOnd==r0.c[i])] = p.adjust(r0$p.val[which(r0$cOnd==r0.c[i])], method = "BH")};rm(i)
    r0$cOnd = NULL
    write.csv(r0, d.nam, quote = F, row.names = F)
    cat("p.adj calculation done:",date(),"\n")
    rm(d.cat,f0,d.Comp, r0)
}

##### Volcano plots #####
d.vol = read.csv(d.nam, header = T)
d.vol = d.vol[which(d.vol$gene %in% GTF$Locus.Tag[GTF$Feature.Type=="CDS"]),]
d.vol$cOnd = paste0(d.vol$cond1,d.vol$cond2)
d.vol.c = unique(d.vol$cOnd)

cat("Start plotting Volcanos:",date(),"\n")
for(i in 1:length(d.vol.c)){
    d.t = d.vol[which(d.vol$cOnd==d.vol.c[i]),]
    d.t = d.t[!is.na(d.t$log2FC),]
    d.t$log2FC[is.infinite(d.t$log2FC)] = ifelse(d.t$log2FC[is.infinite(d.t$log2FC)]>0,max(d.t$log2FC[is.finite(d.t$log2FC)])+1, min(d.t$log2FC[is.finite(d.t$log2FC)])-1)
    cat("Processing scenario",unique(d.t$cond1),"vs",unique(d.t$cond2),"(",i,"/",length(d.vol.c),") :",date(),"          \r")
    pdf(paste0(pT[2],"dNdS_volcano_",unique(d.t$cond1),"_",unique(d.t$cond2),"-vol.pdf"), width = 10, height = 12)
    print(EnhancedVolcano(d.t,
                lab = d.t$gene,
                x = "log2FC",
                y = 'p.val',
                xlab = bquote(~Log[2] ~ "[" ~ Log[10] ~ "median ratio] change"),
                title = paste0(capFirst(unique(d.t$cond1))," (L) vs ",capFirst(unique(d.t$cond2))," (R)"),
                pCutoff= .1, #FDR 
                pCutoffCol = "p.adj", # EnhancedVolcano suggestion: false-discovery rate
                FCcutoff = 1,
                pointSize = 2,
                labSize = 0,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendLabels = c('NS', expression(Log[2] ~ LMRC), 'Adjusted p-value', expression(adjusted ~ p - value ~ and ~ log[2] ~ LMRC)),
                legendPosition = 'bottom',
                legendLabSize = 12,
                legendIconSize = 4,
                drawConnectors = F,
#                max.overlaps = Inf,
                widthConnectors = 0))
    invisible(dev.off())
};rm(i, d.t)
cat("\nDone plotting Volcanos:",date(),"\n")
