#!/bin/bash
# author: ph-u
# script: gFlanking.sh [from PREP_genomeFlanking.sh]
# desc: clip flanking regions of every gene of a pair of cds & genomic FASTA files
# in: bash PREP_genomeFlanking.sh [../relative/2/cds.fa] [../path/2/genomic.fa] [(optional) flanking region length]
# out: [../relative/2/cds]-flanking.csv
# arg: 2 | 3
# date: 20240121, 20241120, 20250223

##### WARNING #####
# Due to difference in cds genome header formats, this program might require manual modification before use
###################

[[ -z $2 ]]&&head -n 5 $0 | tail -n 1&&exit
cDs=$1; gNm=$2
[[ -z $3 ]]&&lEn=100||lEn=$3 # default including upstream xxx sequence (GA-rich region) & downstream termination sequence (T/G-rich, variable length & location from the gene)
pT=`dirname ${cDs}`
bNam=`basename ${cDs} | sed -e "s/[.]fa//"`
[[ -f ${pT}/${bNam}-fT.csv ]]&&rm ${pT}/${bNam}-fT.csv
[[ -f ${pT}/${bNam}-locus.csv ]]&&rm ${pT}/${bNam}-locus.csv

grep -e ">" ${cDs} | sed -e "s/\[locus/@/" | cut -f 2 -d "@" | sed -e "s/=/@/" | cut -f 2 -d "@" | sed -e "s/\]/@/" -e "s/\[location/@/" | cut -f 1,3 -d "@" | sed -e "s/=/@/" -e "s/\]/@/" | cut -f 1,3 -d "@" | sed -e "s/[.][.]/@/" > ${pT}/${bNam}-locus.csv

i0=`wc -l < ${pT}/${bNam}-locus.csv`
sqNm=`basename ${gNm} | sed -e "s/[.]f/@/" | cut -f 1 -d "@"`
mxGen=`tail -n 1 ${gNm} | wc -c`
date
for i in `seq 1 ${i0}`;do
    printf "${i} / ${i0} `date`       \r"
    i2=`head -n ${i} ${pT}/${bNam}-locus.csv | tail -n 1`
    gNam=`echo -e "${i2}" | cut -f 1 -d "@"`
    sqSt=`echo -e "${i2}" | cut -f 2 -d "@" | sed -e "s/complement(//" -e "s/>//g" -e "s/<//g"`
    sqEd=`echo -e "${i2}" | cut -f 3 -d "@" | sed -e "s/)//g" -e "s/>//g" -e "s/<//g"`
    if [[ `echo -e "${sqSt}" | grep -e "join" | wc -l` -gt 0 ]];then
        sqSt=`echo -e "${sqSt}" | sed -e "s/join(//"`
        sqEd=`echo -e "${sqEd}" | rev | sed -e "s/[.][.]/@/" | cut -f 1 -d "@" | rev`
    fi
    if [[ $(( ${sqSt} -1 )) -gt ${lEn} ]];then
        fBef=`head -n 2 ${gNm} | tail -n 1 | cut -c $(( ${sqSt} -${lEn} ))-$(( ${sqSt} -1 ))`
    elif [[ ${sqSt} -le 1 ]];then
        fBef=`head -n 2 ${gNm} | tail -n 1 | cut -c $(( ${mxGen} -${lEn} ))-${mxGen}`
    else
        fBef=`head -n 2 ${gNm} | tail -n 1 | cut -c 1-$(( ${sqSt} -1 ))`
        fLen=`echo -e "${fBef}" | wc -c`
        fB2=`head -n 2 ${gNm} | tail -n 1 | cut -c $(( ${mxGen} -${fLen} ))-${mxGen}`
        fBef=`echo -e "${fB2}${fBef}"`
    fi
    if [[ $(( ${sqEd} +${lEn} )) -lt ${mxGen} ]];then
        fAft=`head -n 2 ${gNm} | tail -n 1 | cut -c $(( ${sqEd} +1 ))-$(( ${sqEd} +${lEn} ))`
    elif [[ ${sqEd} -ge ${mxGen} ]];then
        fAft=`head -n 2 ${gNm} | tail -n 1 | cut -c 1-${lEn}`
    else
        fAft=`head -n 2 ${gNm} | tail -n 1 | cut -c $(( ${sqEd} +1 ))-${mxGen}`
        fLen=`echo -e "${fAft}"`
        fA2=`head -n 2 ${gNm} | tail -n 1 | cut -c 1-$(( ${lEn} -${fLen} ))`
        fAft=`echo -e "${fAft}${fA2}"`
    fi
    [[ ${fBef} == "" ]]&&echo -e "WARNING: cannot find pre-flanking region\nWARNING: Check genome sequence was captured for 'fBef' variable"&&exit
    [[ ${fAft} == "" ]]&&echo -e "WARNING: cannot find post-flanking region\nWARNING: Check genome sequence was captured for 'fAft' variable"&&exit
    echo -e "${sqNm},${gNam},${sqSt},${sqEd},${fBef},${fAft}" >> ${pT}/${bNam}-fT.csv
done
printf "\n`date`\n"

mv ${pT}/${bNam}-fT.csv ${pT}/${bNam}-flanking.csv
rm ${pT}/${bNam}-locus.csv
exit
