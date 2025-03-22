#!/bin/env bash
# author: ph-u
# script: getPA14Seq.sh
# desc: snippet getting PA14 equivalent from PseudoCap seq
# in: bash getPA14Seq.sh 'ATCG'
# out: stdout
# arg: 1
# date: 20250322

i=$1
i1=`echo -e "${i}" | sed -e "s/ //g"`
i0=`grep -n ${i1} ../data/GCF0000146251_cds.fa | cut -f 1 -d ":"`
echo -e "\n`grep -n ${i1} ../data/GCF0000146251_cds.fa | wc -l`"
head -n $(( ${i0} -1 )) ../data/GCF0000146251_cds.fa | tail -n 1
echo -e "\n"
exit
