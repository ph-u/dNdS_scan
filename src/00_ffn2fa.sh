#!/bin/bash
# author: ph-u
# script: 00_ffn2fa.sh
# desc: convert multiline ffn protein-coding annotated genes to FASTA format
# in: bash 00_ffn2fa.sh [../relative/path/2/basename.ffn]
# out: raw/[basename].fa
# arg: 1
# date: 20231220

[[ -z $1 ]]&&head -n 5 $0 | tail -n 1&&exit

i=$1
if [[ `echo -e "${i}" | rev | cut -f 1 -d "." | rev` == "ffn" ]];then
    nAm=`sed -e "s/.ffn//" <<< ${i}`
elif [[ `echo -e "${i}" | rev | cut -f 1 -d "." | rev` == "fna" ]];then
    nAm=`sed -e "s/.fna/_genomic/" <<< ${i}`
fi
echo -e ${nAm}
[[ ${nAm} == "" ]]&&exit
[[ -f ${nAm}_t.txt ]]&& rm ${nAm}_t.txt
[[ -f ${nAm}_dup.txt ]]&&rm ${nAm}_dup.txt
grep -n ">" ${i} | cut -f 1 -d ":" >> ${nAm}_t.txt
echo -e "$(( `wc -l < ${i}` +1 ))" >> ${nAm}_t.txt
i0=`wc -l < ${nAm}_t.txt`
[[ -f ${nAm}.fa ]]&&rm ${nAm}.fa
touch ${nAm}.fa
echo -e "Multi- to single-line seq - `date`"
for i1 in `seq 2 ${i0}`;do
    printf "${i1} / ${i0} - `date`       \r"
    i2=`head -n ${i1} ${nAm}_t.txt | tail -n 1` # seq end line
    i3=`head -n $(( ${i1} -1 )) ${nAm}_t.txt | tail -n 1` # desc line
    head -n $(( ${i2} -1 )) ${i} | tail -n +${i3} > ${nAm}_t0.txt
    i4=`head -n 1 ${nAm}_t0.txt | cut -f 1 -d " "`
    if [[ `grep -e "${i4}" ${nAm}.fa | wc -l` -eq 0 ]];then # exclude curation duplications
        head -n 1 ${nAm}_t0.txt >> ${nAm}.fa
        sEq=`tail -n +2 ${nAm}_t0.txt`
        echo "${sEq//[$'\n']}" >> ${nAm}.fa # http://stackoverflow.com/questions/19345872/ddg#19347380
    else
        head -n 1 ${nAm}_t0.txt >> ${nAm}_dup.txt
    fi
done; printf "\n"
[[ -f ${nAm}_t.txt ]]&& rm ${nAm}_t.txt
[[ -f ${nAm}_t0.txt ]]&& rm ${nAm}_t0.txt
echo -e "Finish process - `date`"
exit
