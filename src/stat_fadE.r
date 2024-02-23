#!/bin/env Rscript
# author: ph-u
# script: stat_fadE.r
# desc: statistical analysis of domain-specific dN/dS values
# in: Rscript stat_fadE.r
# out: res/stat_fadE.csv
# arg: 0
# date: 20240221

pT = paste0("../",c("data","res"),"/")
fadE1 = read.csv(paste0(pT[1],"01_PAO1_107_PA0506.csv"), header = T)
fadE2 = read.csv(paste0(pT[1],"01_PAO1_107_PA0508.csv"), header = T)
# b-sheet-domain: 176-280
# N-domain: 45-162
# CI-domain: 6-45, 94-117, 283-448
# CII-domain: 450-586

dOmain = function(pOx){pOs = 1+floor(pOx/3)
    if(pOs >= 176 & pOs <= 280){return("b")}else if(pOs >= 45 & pOs <= 162){return("n")}else if(pOs >= 450 & pOs <= 586){return("c2")}else if((pOs >= 6 & pOs <= 45) | (pOs >= 94 & pOs <= 117) | (pOs >= 283 & pOs <= 448)){return("c1")}else{return("z")}
}
domainRange = function(sT,eD){
    aa = seq(sT,eD,3)
    for(z in 1:length(aa)){aa[z] = dOmain(as.numeric(aa[z]))}
    a0 = numeric(4)
    if("b" %in% aa){a0[1] = 1}
    if("n" %in% aa){a0[2] = 1}
    if("c1" %in% aa){a0[3] = 1}
    if("c2" %in% aa){a0[4] = 1}
    return(a0)
}

rEs0 = c("B","N","C1","C2")
fadE1.0 = as.data.frame(matrix(NA, nr = nrow(fadE1), nc = length(rEs0)))
fadE2.0 = as.data.frame(matrix(NA, nr = nrow(fadE2), nc = length(rEs0)))
colnames(fadE1.0) = colnames(fadE2.0) = rEs0
i0 = max(c(nrow(fadE1),nrow(fadE2)))

cat(date(),"\n")
for(i in 1:i0){
    cat(i,"/",i0,"-",round(i/i0,2)*100,"% -",date(),"       \r")
    if(i <= nrow(fadE1)){fadE1.0[i,] = domainRange(fadE1$start[i], fadE1$end[i])}
    if(i <= nrow(fadE2)){fadE2.0[i,] = domainRange(fadE2$start[i], fadE2$end[i])}
};rm(i)
cat("\n",date(),"\n")

write.csv(fadE1.0, paste0(pT[1],"stat_fadE1.csv"), row.names = F, quote = F)
write.csv(fadE2.0, paste0(pT[1],"stat_fadE2.csv"), row.names = F, quote = F)
