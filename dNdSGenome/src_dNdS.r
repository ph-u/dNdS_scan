#!/bin/env Rscript
# author: ph-u
# script: src_dNdS.r
# desc: environment for dN/dS calculations
# in: source("src_dNdS.r")
# out: NA
# arg: 0
# date: 20230317, 20230324 (nucleotide ambiguity), 20231217 (implement as a module)

#https://bioinformatics.cvr.ac.uk/calculating-dnds-for-ngs-datasets/
rEf = read.csv("codon.csv", header=T)
nT = read.csv("nt_amb.csv", header=T)
rEf$cDNA = gsub("u","t",rEf$RNA) # DNA2protein translation
nUcleotide = c("a","t","c","g")
cOmplement = read.csv("nt_complementary.csv", header=T)

##### functions #####
nFlip = function(x="a",sTrand="-",x0=cOmplement){return(if(sTrand=="-"){x0[which(x0[,1]==tolower(x)),2]}else{tolower(x)})}
codon2aa = function(cOdon="atg", dIct = rEf){return(dIct$aa[grep(tolower(cOdon),dIct$cDNA)])} ## codon 2 protein

n2p = function(cOdon="atg", n = nUcleotide, nT0 = nT){
	if(all(strsplit(cOdon,"")[[1]] %in% n)){return(codon2aa(cOdon))}
	else{
		n0 = strsplit(cOdon,"")[[1]]
		if(any(n0=="z")){return("-")}; n9 = ";"
		n1 = which(n0 %in% nT$code)
		n2 = 1; for(i in 1:length(n1)){n2 = n2*length(strsplit(nT0$nucleotide[which(nT0$code==n0[n1[i]])],n9)[[1]])}
		n3 = vector(mode="list", length=length(n0)); for(i in 1:length(n3)){if(n0[i] %in% n){n3[[i]] = n0[i]}else{n3[[i]] = strsplit(nT0$nucleotide[which(nT0$code==n0[i])],n9)[[1]]}}
		n4 = as.data.frame(matrix(nr=n2,nc=2))
		n5 = 1;for(i0 in n3[[1]]){for(i1 in n3[[2]]){for(i2 in n3[[3]]){n4[n5,1] = paste0(c(i0,i1,i2), collapse="");n5 = n5+1}}}
		for(i in 1:nrow(n4)){n4[i,2] = codon2aa(n4[i,1])}
		n5 = unique(n4[,2])
		if(length(n5)>1){n5 = paste0("[",paste0(n5[order(n5)],collapse=""),"]")}
		return(n5)
	}}

nt2prot = function(ntSeq="atgggt"){
    x0 = strsplit(tolower(ntSeq), "")[[1]]
    xSeq = seq(1,length(x0),3)
    x1 = numeric(length(xSeq))
    for(i in 1:length(xSeq)){
        x1[i] = n2p(paste0(x0[xSeq[i]],x0[xSeq[i]+1],x0[xSeq[i]+2]))
        if(nchar(x1[i])>1){x1[i] = substr(x1[i],2,2)} # if amino acid residue is ambiguious, select one out (enable visualization, but could introduce protein structure error!!!)
        if(x1[i]==9){break} }
    return(paste0(x1[1:(i-1)], collapse=""))
}

nRp = function(cOdon="atg", n = "c", sIte = 1){ ## codon nt replacement
	cOdon = tolower(cOdon);n = tolower(n)
	if(sIte<2){return(paste0(n,substr(cOdon,2,3)))}
	else if(sIte>2){return(paste0(substr(cOdon,1,2),n))}
	else{return(paste0(substr(cOdon,1,1),n,substr(cOdon,3,3)))}}

nsDet = function(c1="atg",c2="att"){if(nchar(n2p(tolower(c1)))>1 | nchar(n2p(tolower(c2)))>1){return(c(NA,NA))}else if(n2p(tolower(c1))==n2p(tolower(c2))){return(c(0,1))}else{return(c(1,0))}} # c(non-synonymous, synonymous)

nSite = function(cOdon="atg", site=1, nT=nUcleotide){ ## single position missense
	x1 = c(); for(i in nT[nT!=substr(cOdon,site,site)]){
		if(site<2){x0 = paste0(i,substr(cOdon,2,3))}
		else if(site>2){x0 = paste0(substr(cOdon,1,2),i)}
		else{x0 = paste0(substr(cOdon,1,1),i,substr(cOdon,3,3))}
		x1 = c(x1, n2p(x0))}
	return(sum(x1!=n2p(cOdon))/3)}

codonNsite = function(cOdon="atg"){ ## codon missense vector
	x2 = c(); for(i in 1:3){x2 = c(x2, nSite(cOdon,i))}
	return(c(sum(x2),3-sum(x2)))}

seqNsite = function(sEq = "atgaaacccgggttttaaa"){ ## gene sequence missense vector: c(N, S)
	i0 = seq(1,nchar(sEq)-nchar(sEq)%%3,3)
	x3 = as.data.frame(matrix(0,nr=,nc=2))
	for(i in 1:length(i0)){x3[i,] = codonNsite(substr(sEq,i0[i],i0[i]+2))}
	return(as.numeric(colSums(x3)))}

codonMTpath = function(oRi,fIn){
	x0 = strsplit(c(tolower(oRi),tolower(fIn)),"")
	i3 = which(x0[[1]] != x0[[2]]) # identify nt difference(s)
	i0 = length(i3) # count number of nt difference
	if(i0==0){return(c(0,0))}else if(i0==1){nsDet(oRi,fIn)}
	else if(i0==2){ # 2 nt differences, 2 mutation paths
		x1 = c(nRp(oRi,x0[[2]][i3[1]],i3[1]), nRp(oRi,x0[[2]][i3[2]],i3[2]))
		x2 = c(0,0);for(i in 1:length(x1)){x2 = x2 + nsDet(oRi,x1[i]) + nsDet(x1[i],fIn)}
		return(x2/i0)}
	else{x1 = as.data.frame(matrix(nr=factorial(i0),nc=i0-1))
		x1[,1] = rep(c(nRp(oRi,x0[[2]][1],1),nRp(oRi,x0[[2]][2],2),nRp(oRi,x0[[2]][3],3)),each=2)
		x2 = c(2,3,1,3,1,2); x3 = c(0,0)
		for(i in 1:nrow(x1)){x1[i,2] = nRp(x1[i,1],x0[[2]][x2[i]],x2[i])
			x3 = x3 + nsDet(oRi,x1[i,1]) + nsDet(x1[i,1],x1[i,2]) + nsDet(x1[i,2],fIn)}
		return(x3/nrow(x1))}}

samNS = function(sAmple = "atgaaacgcggctactaaa", rEference = "atgaaacccgggttttaaa"){ ## seq compare: c(Nd, Sd)
	i0 = min(nchar(sAmple),nchar(rEference)) # shorter sample as length reference
	i1 = seq(1,i0-i0%%3,3) # codon start positions
	i2 = which(strsplit(tolower(sAmple),"")[[1]] != strsplit(tolower(rEference),"")[[1]]) # single nucleotide polymorphisms (SNP)
	if(length(i2)>0){
		i3 = c();for(i in 1:length(i2)){i3 = c(i3,max(i1[i1<=i2[i]]))};i3 = table(i3) # codon containing SNPs
		x = c("codonPos","mutation","refCodon","mutCodon","reference","mutant")
		x0 = as.data.frame(matrix(nr=length(i3), nc=length(x))) # record protein change
		colnames(x0) = x
		x = c(0,0); for(i in 1:length(i3)){ # N-dist (Nd), S-dist (Sd) calculation
			x1 = c(substr(rEference,names(i3)[i],as.numeric(names(i3)[i])+2), substr(sAmple,names(i3)[i],as.numeric(names(i3)[i])+2))
			x0[i,] = c(names(i3)[i],i3[i],x1, n2p(x1[1]), n2p(x1[2]))
			x = x + codonMTpath(x1[1],x1[2])
		};x0$type = ifelse(x0$reference==x0$mutant,"Synonymous",ifelse(x0$mutant==9,"Nonsense","Missense"))
	}else{x = c(0,0); x0 = "Sequence identical"}
	return(list(x,x0))}

dNdS.rt = function(sAm = "atgaaacgcggctactaaa", rEf = "atgaaacccgggttttaaa"){ ## dN/dS calculation
	NdSd = samNS(sAm,rEf);NS = seqNsite(rEf)
	p = NdSd[[1]]/NS
	dNdS0 = ifelse(any(p>.75) | any(is.na(p)),NA,log(1-4*p[1]/3)/log(1-4*p[2]/3)) #full eq: (-.75*log(1-4*p[1]/3))/(-.75*log(1-4*p[2]/3))
	x = c(dNdS0,p,NdSd[[1]],NS);names(x) = c("dNdS","pN","pS","Nd","Sd","N","S")
	return(list(x,NdSd[[2]]))}
