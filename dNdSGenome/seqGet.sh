#!/bin/bash
# author: ph-u
# script: seqGet.sh
# desc: extract blast matching sequence
# in: bash seqGet.sh [strain_prefix] [gene_name] [database_name]
# out: appending onto data/00_[strain_prefix]_[gene_name]_db.fa
# arg: 3
# date: 20240125

sTrain=$1; gEne=$2; iPcd=$3

grep -v "^$" ../res/${sTrain}_${gEne}.txt > ../res/${sTrain}_${gEne}_t.txt
mv ../res/${sTrain}_${gEne}_t.txt ../res/${sTrain}_${gEne}.txt

if [[ `grep -e "No hits found" ../res/${sTrain}_${gEne}.txt | wc -l` -eq 0  ]];then
##### capture best match result #####
    if [[ `grep -e "Score =" ../res/${sTrain}_${gEne}.txt | wc -l` -lt 2  ]];then
        lN=$(( `grep -n "Lambda" ../res/${sTrain}_${gEne}.txt | head -n 1 | cut -f 1 -d ":"` -1  ))
    else
        lN=$(( `grep -n "Score =" ../res/${sTrain}_${gEne}.txt | head -n 2 | tail -n 1 | cut -f 1 -d ":"` -1  ))
    fi
    lN0=`grep -n "Score =" ../res/${sTrain}_${gEne}.txt | head -n 1 | cut -f 1 -d ":"`
    head -n ${lN} ../res/${sTrain}_${gEne}.txt | tail -n +${lN0} > ../data/00_${sTrain}_${gEne}_t.txt # capture the best match gene content

##### get best match database sequence #####
    grep -e "Sbjct" ../data/00_${sTrain}_${gEne}_t.txt | tr -s " " > ../data/00_${sTrain}_${gEne}_t0.txt
    dbSt=`head -n 1 ../data/00_${sTrain}_${gEne}_t0.txt | tail -n 1 | cut -f 2 -d " "`
    dbEd=`tail -n 1 ../data/00_${sTrain}_${gEne}_t0.txt | cut -f 4 -d " "`
    [[ ${dbSt} -gt ${dbEd}  ]]&&pM="minus"||pM="plus"

    echo -e ">${iPcd};${dbSt}-${dbEd};${pM}" >> ../data/00_${sTrain}_${gEne}_db.fa

    dbSeq=`cat ../data/00_${sTrain}_${gEne}_t0.txt | cut -f 3 -d " "`
    echo "${dbSeq//[$'\n']}" >> ../data/00_${sTrain}_${gEne}_db.fa
    [[ -f ../data/00_${sTrain}_${gEne}_t0.txt ]] && rm ../data/00_${sTrain}_${gEne}_t0.txt
    [[ -f ../data/00_${sTrain}_${gEne}_t.txt ]] && rm ../data/00_${sTrain}_${gEne}_t.txt
fi
exit