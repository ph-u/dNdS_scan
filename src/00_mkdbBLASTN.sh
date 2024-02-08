#!/bin/bash
# author: ph-u
# script: 00_mkdbBLASTN.sh
# desc: section sequence for blastn local query
# in: bash mkdbBLASTN.sh [../relative/path/2/sourceFasta] [../relative/path/2/dbDirectory]
# out: [../relative/path/2/dbDirectory]*
# arg: 2
# date: 20230328, 20231223

##### ref #####
# https://ncbi.github.io/magicblast/cook/blastdb.html
[[ -z $2 ]]&&head -n 5 $0 | tail -n 1&&exit

##### env #####
i0=`pwd`
i1="$HOME/Desktop/ncbi_cli/ncbi-blast-2.13.0+/bin"; i2=$2
i4=$1
echo -e "$PATH" > tmp.txt
[[ `grep -e "${i1}" tmp.txt | wc -l` != 1 ]]&&export PATH="${i1}:$PATH"
rm tmp.txt

##### make blastdb #####
for i in `ls ${i4}/*.fna`;do
    i5=`echo -e "${i}" | rev | cut -f 1 -d "/" | rev`
    echo -e "${i5}"
    makeblastdb -in ${i} -dbtype nucl -parse_seqids -out ${i2}/${i5} -title ${i5} 1> p_mkdbBLASTNmsg.log
done
exit
