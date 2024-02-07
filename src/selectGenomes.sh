#!/bin/bash
# author: ph-u
# script: selectGenomes.sh
# desc: select E coli genomes from all NCBI genomes
# in: bash selectGenomes.sh [genomes_selected] [(optional) genomes_pool.xxx] [(optional) exclusion_list.xxx]
# out: raw/selectGenomes.csv
# arg: 1 - 3
# date: 20240123

# https://unix.stackexchange.com/questions/506207/fast-way-to-extract-lines-from-a-large-file-based-on-line-numbers-stored-in-anot/506236#506236
[[ -z $1 ]] && head -n 5 $0 | tail -n 1 && exit
[[ -z $2 ]] && gPool="../raw/eColiMetaData.tsv" || gPool=$2
[[ -z $3 ]] && xList="../raw/accession.csv" || xList=$3

echo -e "Counting lines in genome pool - `date`"
i1=`wc -l < ${gPool}`
head -n 1 ${gPool} > ../raw/selectGenomes.csv
i=0
date
while [[ ${i} -lt $1 ]];do
    i0=$(( ${RANDOM} % ${i1} ))
    printf "${i0}     \r"
    if [[ ${i0} -gt 1 ]] && [[ `grep -e ${i0} ${xList} | wc -l` -eq 0 ]];then
        head -n ${i0} ${gPool} | tail -n 1 >> ../raw/selectGenomes.csv
        i=$(( ${i} +1 ))
    fi
done
date
exit
