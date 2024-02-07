#!/bin/bash 
# author: ph-u
# script: seqCollect_c.sh
# desc: child script calling blastn sequence alignment
# in: bash blastn_cRE.sh [$SLURM_ARRAY_TASK_ID]
# out: res/[gene]-[clinical].txt
# arg: 1
# date: 20240125

tID=$1      
i2=`head -n ${tID} iDx.csv | tail -n 1 | cut -f 1 -d ","` # source FASTA file name
i3=`echo -e "${i2}" | rev | cut -f 1 -d "_" | rev` # gene name
i5=`echo -e "${i2}" | rev | sed -e "s/_/@/" | rev | cut -f 1 -d "@"` # strain name

if [[ ! -f ../data/${i2}.fa ]];then
    k0=`grep -e ${i3} ../data/${i5}_meta.txt | cut -f 1 -d ":"`
    head -n $(( ${k0} +1  )) ../data/${i5}.fa | tail -n 2 > ../data/${i2}.fa
fi

tail -n +2 ../data/00_${i2}.csv | grep -v "blastnERR%" | uniq | cut -f 2 -d "," > ../data/00_${i2}-dbList.txt

p0=`wc -l < ../data/00_${i2}-dbList.txt`
[[ -f ../data/00_${i2}_db.fa ]] && rm ../data/00_${i2}_db.fa
touch ../data/00_${i2}_db.fa
for i0 in `seq 1 ${p0}`;do
    i1=`head -n ${i0} ../data/00_${i2}-dbList.txt | tail -n 1` # database
    blastn -query ../data/${i2}.fa -db ${i1} -out ../res/${i5}_${i3}.txt
    bash seqGet.sh ${i5} ${i3} ${i1}
    echo -e "Done ${i1} matched seq extraction - `date`"
done
[[ -f ../data/00_${i2}-dbList.txt ]] && rm ../data/00_${i2}-dbList.txt
exit
