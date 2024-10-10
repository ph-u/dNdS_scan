#!/bin/env Rscript
# author: ph-u
# script: dNdS_median-PCA.r
# desc: clustering dNdS medians by source
# in: Rscript dNdS_median-PCA.r
# out: {data,res}/dNdS_median-PCA.*
# arg: 0
# date: 20240620

source("p_src.r")
source("p_metabolism_PAO1.r")
dNdS.md = paste0(pT[1],"dNdS_median-PCA.csv")
med0.c = unique(mEta$sOurce)

##### Collect median dN/dS by source and gene if file not already existed #####
if(!file.exists(dNdS.md) | file.info(dNdS.md)$size<=126037){
    f0 = list.files(pT[3], "dbSum")
    f = data.frame(fNam = f0, gEne = read.table(text = gsub("_PA","@PA",gsub("--","@",f0)), sep = "@")[,2])
    dNdS.med = as.data.frame(matrix(nr = nrow(f), nc = length(med0.c)*3))
    row.names(dNdS.med) = f$gEne
    colnames(dNdS.med) = paste0(med0.c,rep(c("",".inf",".NA"), each = length(med0.c)))
    dNdS.cat = vector(mode = "list", length = length(med0.c))
    names(dNdS.cat) = med0.c
    for(i in 1:length(dNdS.cat)){dNdS.cat[[i]] = mEta$assemblyInfo.genbankAssmAccession[which(mEta$sOurce==names(dNdS.cat)[i])]};rm(i)
    cat("Starting dNdS median value data-collection:",date(),"\n")
    i0=i1=1;repeat{ cat("Processing gene",row.names(dNdS.med)[i0],"(",i0,"/",nrow(dNdS.med),")", ", source",colnames(dNdS.med)[i1],"(",i1,"/",ncol(dNdS.med),")", date(), "                                       \r")
        if(i1==1){
            d0 = read.csv(paste0(pT[3], f$fNam[which(f$gEne==row.names(dNdS.med)[i0])]), header = T)
            d0$dNdS[is.infinite(d0$dNdS)] = Inf # -Inf in data actually denote non-synonymous only, which should be high positives
            d0$clinical = read.table(text = gsub("_ASM","@",d0$clinical), sep = "@")[,1]
        }
        med1 = d0$dNdS[which(d0$clinical %in% dNdS.cat[[colnames(dNdS.med)[i1]]])]
        dNdS.med[i0,i1+length(med0.c)*0:2] = c(median(med1, na.rm = T),sum(is.infinite(med1)), sum(is.na(med1)))
        i1 = i1+1
        if(i1 > length(med0.c)){i1=1;i0=i0+1}
        if(i0 > nrow(dNdS.med)){break}
    };rm(i0,i1,med1);cat("\nDone:",date(),"\n")
    rm(dNdS.cat, f0)
    write.csv(dNdS.med, dNdS.md, row.names = T, quote = F)
};rm(f)

##### dNdS data curation #####
dNdS.med = read.csv(dNdS.md, header = T, row.names = 1)
dNdS.breakpt = rowSums(dNdS.med[,-c(1:length(med0.c))])
dNdS.medBP = is.finite(rowSums(dNdS.med[,c(1:length(med0.c))]))
dNdS.med = dNdS.med[which(dNdS.medBP & dNdS.breakpt<(nrow(mEta)*.6)),1:length(med0.c)] # at least 40% data with valid numbers; 0.7 -> 5059 genes; 0.6 -> 4792 genes (.5 / .6 optimal; small number = over-fitting)
dMax = max(dNdS.med, na.rm = T)+1
for(i in 1:ncol(dNdS.med)){dNdS.med[is.infinite(dNdS.med[,i]),i] = dMax};rm(i)
#dNdS.med = dNdS.med[which(row.names(dNdS.med) %in% gBioc$protList),]

##### Colour grouping #####
cat("Categorizing colour groups:",date(),"\n")
#dNdS.med$gP = gBioc$metabolism[match(row.names(dNdS.med), gBioc$protList)]

## PAO1 features
#fEat = unique(read.csv(paste0(pT[4], "features.txt"), header = T, sep = "\t", comment.char = "", quote = "")[,c("Locus.Tag", "Gene.Length","Estimated.MW..kDa.", "Estimated.Isoelectric.Point..pI.", "Estimated.Charge..ph7.")])

## GO.Term
cogLim = c(3,99,200)
cogData = read.table("../raw/gene_ontology_tab.txt", sep = "\t", header = T, comment.char = "", quote = "")
cogData$GO.Term = gsub("[[]",";",gsub("[]]",";",cogData$GO.Term))
cTable = table(cogData$GO.Term)
c0.raw = data.frame(Category = ifelse(cTable>cogLim[1],ifelse(cTable<cogLim[2],names(cTable),ifelse(cTable>cogLim[3],"Ext.Func.Redun","High.Func.Redun")),ifelse(cTable==1,"Uniq.Func.Redun","Low.Func.Redun")), sRc = names(cTable))
row.names(c0.raw) = 1:nrow(c0.raw)

pT01 = "../../1_04_proteomics/raw/Pinyu/"
c0.PY = unique(read.csv(paste0(pT01,"go_filtered_output.csv"), header = T)[,c("Category", "GO.Term")])
c0.PY$GO.Term = gsub("[[]",";",gsub("[]]",";",c0.PY$GO.Term))
c0.PY$Category[grep("/",c0.PY$Category)] = "Others"

for(i in which(!(1:nrow(c0.raw) %in% grep("Func.Redun",c0.raw$Category)))){
    i0 = grep(c0.raw$Category[i], c0.PY$GO.Term)
    c0.raw$Category[i] = ifelse(length(i0)==0, "Others", c0.PY$Category[i0])
};rm(i,i0)

## Secretome vs Matrixome
#c0.mtx = unique(read.csv(paste0(pT01,"go_filtered_matrixome.csv"), header = T)[,c("PANum","Category", "GO.Term")])
#c0.stm = unique(read.csv(paste0(pT01,"go_filtered_secretome.csv"), header = T)[,c("PANum","Category", "GO.Term")])
#c0.raw = unique(data.frame(Category=c(rep("matrixome", length(c0.mtx$PANum)), rep("secretome", length(c0.stm$PANum))), sRc = c(c0.mtx$PANum, c0.stm$PANum)))

#gbCol = data.frame(gP = c(unique(c0.raw$Category),"Multiple","Unknown"), cOl = cBp[1:(2+length(unique(c0.raw$Category)))])

for(i in 1:nrow(dNdS.med)){ cat(i,"/",nrow(dNdS.med),"     \r")
#    i0 = unique(c0.raw$Category[c0.raw$sRc %in% row.names(dNdS.med)[i]])
    i0 = cogData$GO.Term[which(cogData$Locus.Tag==row.names(dNdS.med)[i])]
    i0 = unique(c0.raw$Category[c0.raw$sRc %in% i0])
    dNdS.med$gP[i] = ifelse(length(i0)>1, "Multiple", i0) #gBioc$metabolism[match(row.names(dNdS.med), gBioc$protList)]
};rm(i, i0);cat("\nDone:",date(),"\n")

## Merge dN/dS data with physical-chemical properties
#dNdS.med$Locus.Tag = row.names(dNdS.med)
#dNdS.med = merge(dNdS.med,fEat, by = "Locus.Tag", all.x = T)
#row.names(dNdS.med) = dNdS.med$Locus.Tag
#dNdS.med$Locus.Tag = NULL
#fEat.0 = fEat[which(fEat$Locus.Tag %in% row.names(dNdS.med)),]

## Grouping colours & Rare group determination
#gbCol = unique(gBioc[,-1])
gbCol = unique(dNdS.med$gP)
gbCol = data.frame(Category = gbCol[order(gbCol)], cOl = cBp[1:length(gbCol)])
gbCol[nrow(gbCol)+1,] = c("Unknown", cBp[nrow(gbCol)+1])
gbCol[nrow(gbCol)+1,] = c("Rare Groups", cBp[nrow(gbCol)+1])
#colnames(gbCol) = c("Category", "cOl")

dNdS.med$gP[is.na(dNdS.med$gP)] = "Unknown"
i0 = c();for(i in 1:nrow(dNdS.med)){if(all(is.na(dNdS.med[i,-ncol(dNdS.med)]))){i0 = c(i0,i)}};rm(i)
if(length(i0)>0){dNdS.med = dNdS.med[-i0,]};rm(i0)
dNdS.gpTable = table(dNdS.med$gP)
#dNdS.med$gP[which(dNdS.med$gP %in% names(dNdS.gpTable[which(dNdS.gpTable<3)]))] = gbCol$Category[grep("Rare Groups", gbCol$Category)]
dNdS.med$gP = capFirst(dNdS.med$gP)
row.names(dNdS.med) = paste0(substr(dNdS.med$gP,1,1),".",row.names(dNdS.med))

##### PCA with NA #####
p0 = pcaSwap(dNdS.med[,-ncol(dNdS.med)]) # require no inf
p0.r = pcaMethods::pca(dNdS.med[,-ncol(dNdS.med)], method = "ppca") # require no inf
#p0 = prcomp(dNdS.med, center = T)

g0 = ggbiplot(p0, var.scale=1, groups=dNdS.med$gP, ellipse=T, ellipse.prob=.95, labels = row.names(dNdS.med), labels.size=2.5, varname.size = 3, varname.adjust = c(11,5,8,14,6)) +
#g0 = ggbiplot(p0, var.scale=1, ellipse=T, ellipse.prob=.95, labels = row.names(dNdS.med), labels.size=1, varname.size = 4) +
    scale_color_manual(values=setNames(gbCol$cOl, gbCol$gP)) +
	guides(color=guide_legend(title="Gene Functions")) +
	scale_y_continuous(breaks = seq(-40, 20, 10), limits = c(-40, 20))+
	coord_cartesian(xlim=c(-5,25))+
	theme(legend.position = 'bottom', legend.direction = "vertical",
        	panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
	        legend.title = element_text(size=18),
	        axis.text=element_text(size=16),axis.title=element_text(size=14),
	        plot.margin=margin(t=0,r=0,b=0,l=1)) + theme_bw()
#ggsave(paste0(pT[2],"dNdS_median-PCA.pdf"), plot=pcaLAB(g0,round(p0.r@R2*100,1)), width=7, height=7)
ggsave(paste0(pT[2],"dNdS_median-PCA.pdf"), plot=pcaLAB(g0,round(p0.r@R2*100,1)), width=8, height=5)

## PAO1.bioc[grep("PA1984",PAO1.bioc$Genes),]
## cogData[grep("PA1984",cogData$Locus.Tag),]
