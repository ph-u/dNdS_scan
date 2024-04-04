#!/bin/bash
# author: ph-u
# script: seqGet.sh
# desc: extract blast matching sequence
# in: bash seqGet.sh [strain_prefix] [gene_name] [database_name]
# out: appending onto data/00_[strain_prefix]_[gene_name]_db.fa
# arg: 3
# date: 20240125 20240329

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
    cOntig=`grep -e ">" ../res/${sTrain}_${gEne}.txt | head -n 1 | cut -f 1 -d " " | sed -e "s/>//"`

##### get best match database sequence #####
    grep -e "Sbjct" ../data/00_${sTrain}_${gEne}_t.txt | tr -s " " > ../data/00_${sTrain}_${gEne}_t0.txt
    dbSt=`head -n 1 ../data/00_${sTrain}_${gEne}_t0.txt | tail -n 1 | cut -f 2 -d " "`
    dbEd=`tail -n 1 ../data/00_${sTrain}_${gEne}_t0.txt | cut -f 4 -d " "`
    [[ ${dbSt} -gt ${dbEd}  ]]&&pM="minus"||pM="plus"

    if [[ `grep -e "${iPcd}" ../data/00_${sTrain}_${gEne}_db.fa | wc -l` -eq 0 ]];then
        echo -e ">${iPcd};${cOntig};${dbSt}-${dbEd};${pM}" >> ../data/00_${sTrain}_${gEne}_db.fa
        dbSeq=`cat ../data/00_${sTrain}_${gEne}_t0.txt | cut -f 3 -d " "`
        echo "${dbSeq//[$'\n']}" | sed -e "s/-//g" >> ../data/00_${sTrain}_${gEne}_db.fa
    else
        cHead=$(( `grep -n ">${iPcd};${cOntig};${dbSt}-${dbEd};${pM}" ../data/00_${sTrain}_${gEne}_db.fa | cut -f 1 -d ":"` ))
        head -n ${cHead} ../data/00_${sTrain}_${gEne}_db.fa > ../data/00_${sTrain}_${gEne}_db.fat
        dbSeq=`cat ../data/00_${sTrain}_${gEne}_t0.txt | cut -f 3 -d " "`
        echo "${dbSeq//[$'\n']}" | sed -e "s/-//g" >> ../data/00_${sTrain}_${gEne}_db.fat
        tail -n +$(( ${cHead}+2 )) ../data/00_${sTrain}_${gEne}_db.fa >> ../data/00_${sTrain}_${gEne}_db.fat
        mv ../data/00_${sTrain}_${gEne}_db.fat ../data/00_${sTrain}_${gEne}_db.fa
    fi

    [[ -f ../data/00_${sTrain}_${gEne}_t0.txt ]] && rm ../data/00_${sTrain}_${gEne}_t0.txt
    [[ -f ../data/00_${sTrain}_${gEne}_t.txt ]] && rm ../data/00_${sTrain}_${gEne}_t.txt
else
    echo -e ">${iPcd};noHit\n" >> ../data/00_${sTrain}_${gEne}_db.fa
fi
exit
