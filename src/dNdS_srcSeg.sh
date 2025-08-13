#!/bin/bash
# author: ph-u
# script: dNdS_srcSeg.sh
# desc: Segregate & plot dN/dS values by sources
# in: bash dNdS_srcSeg.sh [gene] [del]
# out: data/[strain]_[gene]_[src]--reCon.csv, res/[gene]_reCon.pdf
# arg: 1
# date: 20240424

# i0=`wc -l < ../data/PAO1_proteins.txt`;for i in `seq 1 ${i0}`;do bash dNdS_srcSeg.sh `head -n ${i} ../data/PAO1_proteins.txt | tail -n 1` 1>> testProt.log;done
# for i in PA0342 PA1307 PA1342 PA1442 PA1553 PA2024 PA2120 PA2166 PA2185 PA2622 PA3565 PA3567 PA3653 PA3674 PA3692 PA3793 PA4459 PA5108 PA5221 PA5259 PA5553;do bash dNdS_srcSeg.sh ${i} 1>> testProt.log;done

[[ -z $1 ]]&&echo -e "\n`head -n 5 $0 | tail -n 2`\n"&&exit

echo -e "Processing: $1 - `date`"
mkdir -p ../res/reConBIN
Rscript dNdS_convert.r $1
Rscript residueDNDS.r $1
[[ -z $2 ]]||rm ../data/*_$1_*--reCon.csv
mv ../res/$1_reCon.pdf ../res/reConBIN/$1_reCon.pdf
exit
