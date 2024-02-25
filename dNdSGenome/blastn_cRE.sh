#!/bin/bash
# author: ph-u
# script: blastn_cRE.sh
# desc: child script calling blastn sequence alignment
# in: bash blastn_cRE.sh [$SLURM_ARRAY_TASK_ID]
# out: res/[gene]-[clinical].txt
# arg: 1
# date: 20231216, 20240107

tID=$1
i2=`head -n ${tID} iDx.csv | tail -n 1 | cut -f 1 -d ","` # source FASTA file name
i3=`echo -e "${i2}" | rev | cut -f 1 -d "_" | rev` # gene name
i5=`echo -e "${i2}" | rev | sed -e "s/_/@/" | rev | cut -f 1 -d "@"` # strain name

##### Determine which run requires a re-run #####
lastSlurm=xxxxxxx # replace with last slurm run job id (auto)
dbNum=`wc -l < dbList.txt`
[[ -f ../data/${i5}_${i3}_iDx.csv ]]&&rm ../data/${i5}_${i3}_iDx.csv
touch ../data/${i5}_${i3}_iDx.csv
grep -n " / ${dbNum}" slurm-${lastSlurm}_${tID}.out > ../data/${i5}_${i3}_t.csv
p0=`wc -l < ../data/${i5}_${i3}_t.csv`

## check each run for success or not
for i in `seq 1 ${dbNum}` ;do
    if [[ `grep -e "PASSED ${i} / ${dbNum}" slurm-${lastSlurm}_${tID}.out | wc -l` -gt 0 ]];then
        echo -e "PASSED ${i} / ${dbNum}"
    else
        if [[ ${i} -lt ${p0} ]];then
            pStr=`grep -e " ${i} / ${dbNum}" ../data/${i5}_${i3}_t.csv | tail -n 1 | cut -f 1 -d ":"`
            i0=`head -n $(( ${pStr} +1  )) slurm-${lastSlurm}_${tID}.out | tail -n 1 | grep -e "Calculating dNdS for:" | wc -l`
            pEnd=$(( ${pStr} +3 +${i0} ))
        elif [[ ${i} -eq ${p0} ]];then
            pStr=`grep -e " ${i} / ${dbNum}" ../data/${i5}_${i3}_t.csv | tail -n 1 | cut -f 1 -d ":"`
            pEnd=$(( `wc -l < slurm-${lastSlurm}_${tID}.out` +1 ))
        else
            pStr=0
        fi # if there is a need to section the run
        if [[ ${pStr} -eq 0 ]];then
            echo -e "${i}" >> ../data/${i5}_${i3}_iDx.csv
        else
            tail -n +${pStr} slurm-${lastSlurm}_${tID}.out | head -n $(( ${pEnd}-${pStr} )) > ../data/${i5}_${i3}_t0.csv
            if [[ `wc -l < ../data/${i5}_${i3}_t0.csv` -eq $(( 3+`grep -e "Calculating dNdS for:" ../data/${i5}_${i3}_t0.csv | wc -l` )) ]] && [[ `tail -n 2 ../data/${i5}_${i3}_t0.csv | head -n 1 | grep -e "Exporting: " | wc -l` -gt 0 ]];then
                echo -e "PASSED ${i} / ${dbNum}"
            else
                echo -e "${i}" >> ../data/${i5}_${i3}_iDx.csv
            fi # if the line of output equals to expected
            rm ../data/${i5}_${i3}_t0.csv
        fi # if calculation has done before
    fi # if last slurm-file check already passed
done
rm ../data/${i5}_${i3}_t.csv

##### Resume run #####
p0=`wc -l < ../data/${i5}_${i3}_iDx.csv`
if [[ ${p0} -gt 0 ]];then
    k0=`grep -e "=${i3};" ../data/${i5}_meta.txt | cut -f 1 -d ":"`
    head -n $(( ${k0} +1 )) ../data/${i5}.fa | tail -n 2 > ../data/${i2}.fa
    for i0 in `seq 1 ${p0}`;do
        i4=`head -n ${i0} ../data/${i5}_${i3}_iDx.csv | tail -n 1` # which run needs to repeat
        i1=`head -n ${i4} dbList.txt | tail -n 1` # which database corresponding to the run
        echo -e "`date` - ${i1}: ${i4} / ${dbNum}"
## remove previous error run content
        sed "/${i1}/d" ../data/00_${i5}_${i3}.csv > ../data/00_${i5}_${i3}_00t.csv
        mv ../data/00_${i5}_${i3}_00t.csv ../data/00_${i5}_${i3}.csv
        sed "/${i1}/d" ../data/01_${i5}_${i3}.csv > ../data/01_${i5}_${i3}_01t.csv
        mv ../data/01_${i5}_${i3}_01t.csv ../data/01_${i5}_${i3}.csv
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
fi
rm ../data/${i5}_${i3}_iDx.csv
exit
