#!/bin/bash
# author: ph-u
# script: 00_blastn.sh
# desc: SLURM-optimized blastn alignment pipeline
# in: bash 00_blastn.sh [../relative/path/2/overall_genes.fa]
# out: NA
# arg: 1
# date: 20231216

##### Warning fasta format difference #####
# ln30: ; -> ]
###########################################

[[ -z $1 ]]&&echo&&head -n 5 $0 | tail -n 1&&echo&&exit
rAw=$1
bNam=`echo -e "${rAw}" | rev | cut -f 1 -d "/" | rev | sed -e "s/[.]fa//g"`
mkdir -p ../{data,res}

##### Segregate fasta to every gene #####
printf "Segregating genes from ${rAw} (`date`)"
[[ -f ../data/${bNam}_*.fa ]]&&rm ../data/${bNam}_*.fa
for i in `grep -n ">" ${rAw} | cut -f 1 -d ":"`;do
    if [[ `head -n ${i} ${rAw} | tail -n 1 | grep -e "locus_tag" | wc -l` -gt 0 ]];then
        i0=`head -n ${i} ${rAw} | tail -n 1 | sed -e "s/locus_tag=/@/" | cut -f 2 -d "@" | sed -e "s/]/@/" | cut -f 1 -d "@"`
    else
        i0=`head -n ${i} ${rAw} | tail -n 1 | sed -e "s/> //" -e "s/>//" -e "s/[|]/-/g" | cut -f 1 -d " "`
    fi
    head -n $(( ${i}+1 )) ${rAw} | tail -n 2 > ../data/${bNam}_${i0}.fa
done
printf " -- Done (`date`)\n"

##### list pipeline subprocess details #####
printf "Listing subprocess (`date`)"
[[ -f ${bNam}_iDx.csv ]]&&rm ${bNam}_iDx.csv
for i in `ls ../data/${bNam}_*.fa | grep -v "genom"`;do
    i0=`echo -e ${i} | rev | cut -f 1 -d "/" | rev | sed -e "s/[.]fa//g"`
    echo -e "${i0}," >> ${bNam}_iDx.csv
done
printf " -- Done (`date`)\n"

##### Form child script & execute #####
printf "Assembling run-script 'hpc_c.sh' (`date`)"
cat hpcHead.sh > ${bNam}_hpc_c.sh
p0=`echo -e "$ SLURM_ARRAY_TASK_ID" | sed -e "s/ //"`
echo -e "#SBATCH -J ${bNam}\n#SBATCH --array=1-`wc -l < ${bNam}_iDx.csv`\n\nbash ${bNam}_blastn_c.sh ${p0}" >> ${bNam}_hpc_c.sh
printf " -- Done (`date`)\n"

printf "Modifying child-script 'blastn_c.sh' to ${bNam}_blastn_c.sh (`date`)"
sed -e "s/iDx/${bNam}_iDx/g" blastn_c.sh > ${bNam}_blastn_c.sh
printf " -- Done (`date`)\n"
sbatch ${bNam}_hpc_c.sh
exit
