#!/bin/bash
# author: ph-u
# script: PREP_genomeFlanking.sh
# desc: clip flanking regions of every gene of a pair of cds & genomic FASTA files
# in: bash PREP_genomeFlanking.sh [../relative/2/cds.fa] [../path/2/genomic.fa] [(optional) flanking region length]
# out: [../relative/2/cds]-flanking.csv
# arg: 2 | 3
# date: 20240121, 20241120

##### WARNING #####
# Due to difference in cds genome header formats, this program might require manual modification before use
###################

[[ -z $2 ]]&&head -n 5 $0 | tail -n 1&&exit
cDs=$1; gNm=$2
[[ -z $3 ]]&&lEn=100||lEn=$3 # default including upstream xxx sequence (GA-rich region) & downstream termination sequence (T/G-rich, variable length & location from the gene)
pT=`dirname ${cDs}`
bNam=`basename ${cDs} | sed -e "s/[.]fa//"`
[[ -f ${pT}/${bNam}-fT.csv ]]&&rm ${pT}/${bNam}-fT.csv

i0=`wc -l < ${bNam}_iDx.csv`
sqNm=`basename ${gNm} | sed -e "s/[.]f/@/" | cut -f 1 -d "@"`
date
for i in `seq 1 ${i0}`;do
    printf "${i} / ${i0} `date`       \r"
    gNam=`head -n ${i} ${bNam}_iDx.csv | tail -n 1 | cut -f 1 -d ","` # | rev | cut -f 1 -d "_" | rev`
    i2=`grep -e ${gNam} ${cDs} | sed -e "s/location=/@/" | cut -f 2 -d "@" | sed -e "s/[]]/@/" | cut -f 1 -d "@" | sed -e "s/[()]/@/g" | cut -f 2 -d "@" | sed -e "s/[.][.]/@/"`
    sqSt=`echo -e "${i2}" | cut -f 1 -d "@"`
    sqEd=`echo -e "${i2}" | cut -f 2 -d "@"`
    fBef=`head -n 2 ${gNm} | tail -n 1 | cut -c $(( ${sqSt} -${lEn}  ))-$(( ${sqSt} -1  ))`
    fAft=`head -n 2 ${gNm} | tail -n 1 | cut -c $(( ${sqEd} +1 ))-$(( ${sqEd} +${lEn} ))`
    [[ ${fBef} == "" ]]&&echo -e "WARNING: cannot find pre-flanking region\nWARNING: Check genome sequence was captured for 'fBef' variable"&&exit
    [[ ${fAft} == "" ]]&&echo -e "WARNING: cannot find post-flanking region\nWARNING: Check genome sequence was captured for 'fAft' variable"&&exit
    echo -e "${sqNm},${gNam},${sqSt},${sqEd},${fBef},${fAft}" >> ${pT}/${bNam}-fT.csv
done
printf "\n`date`\n"

mv ${pT}/${bNam}-fT.csv ${pT}/${bNam}-flanking.csv
exit
