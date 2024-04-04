#!/bin/bash
# author: ph-u
# script: PAO1_proteinList.sh
# desc: extract protein collection of PAO1 from literature summary
# in: bash PAO1_proteinList.sh
# out: data/PAO1_proteins.txt
# arg: 0
# date: 20240316

x=`cat ../p_pdf/PAO1_proteomics/PAO1_proteins.csv | sed -e "s/;/,/g"`
x=${x//,/$'\n'}
echo -e "${x}" | grep -v "^$" | grep -e "^PA" | sort -u > ../data/PAO1_proteins.txt
exit
