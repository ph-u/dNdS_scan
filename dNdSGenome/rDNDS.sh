#!/bin/bash
# author: ph-u
# script: rDNDS.sh
# desc: child script to call dN/dS calculation script
# in: bash rDNDS.sh [$SLURM_ARRAY_TASK_ID]
# out: NA
# arg: 1
# date: 20240119

tID=$1
srcFile=`head -n ${tID} iDx.csv | tail -n 1`
pF=`basename ${srcFile} | rev | sed -e "s/_/@/" | rev | cut -f 1 -d "@" | sed -e "s/01_//"` # strain name
gN=`basename ${srcFile} | rev | cut -f 1 -d "_" | rev | sed -e "s/[.]csv//"` # gene name
lC=`grep -e ${gN} ../data/${pF}-flanking.csv | cut -f 2 -d ","` # locus tag

Rscript pairwiseSimilarity.r ${tID}
Rscript rDNDS.r ${srcFile} ${pF} ${gN} ${lC}
exit
