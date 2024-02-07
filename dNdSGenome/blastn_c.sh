#!/bin/bash
# author: ph-u
# script: blastn_c.sh
# desc: child script calling blastn sequence alignment
# in: bash blastn_c.sh [$SLURM_ARRAY_TASK_ID]
# out: res/[gene]-[clinical].txt
# arg: 1
# date: 20231216

tID=$1
i2=`head -n ${tID} iDx.csv | tail -n 1 | cut -f 1 -d ","` # source FASTA file name
i3=`echo -e "${i2}" | rev | cut -f 1 -d "_" | rev` # gene name
i5=`echo -e "${i2}" | rev | sed -e "s/_/@/" -e "s/[|]/@/g" | rev | cut -f 1 -d "@"` # strain name
[[ -f ../data/00_${i5}_${i3}.csv ]]&&rm ../data/00_${i5}_${i3}.csv
[[ -f ../data/01_${i5}_${i3}.csv ]]&&rm ../data/01_${i5}_${i3}.csv
echo -e "gene,clinical,locus,database_start,sample_start,len,flank_befPc,flank_aftPc,dbSeqDir" >> ../data/00_${i5}_${i3}.csv
echo -e "gene,clinical,locus,start,end,flank_befPc,flank_aftPc,dNdS,pN,pS,Nd,Sd,N,S" >> ../data/01_${i5}_${i3}.csv
dbNum=`wc -l < dbList.txt`
for i4 in `seq 1 ${dbNum}`;do
    i1=`head -n ${i4} dbList.txt | tail -n 1`
    echo -e "`date` - ${i1}: ${i4} / ${dbNum}"
## blast
    blastn -query ../data/${i2}.fa -db ${i1} -out ../res/${i5}_${i3}.txt
    bash seqGet.sh ${i5} ${i3} ${i1}
    bash sumBLASTN.sh ${i5} ${i3} ${i4} >> ../data/00_${i5}_${i3}.csv # input geneID, dbID
## rolling dNdS
    Rscript rollingdNdS.r ${i5} ${i3}
    tail -n +2 ../data/01_${i5}_${i3}_t.csv >> ../data/01_${i5}_${i3}.csv
done
rm ../data/${i2}.fa
rm ../res/${i5}_${i3}.txt
rm ../data/01_${i5}_${i3}_t.csv
exit
