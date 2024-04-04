#!/bin/bash
# author: ph-u
# script: 01_dbSum.sh
# desc: Summarize blastn collection on one gene
# in: bash 01_dbSum.sh ../data/[basename].fa
# out: NA
# arg: 1
# date: 20240331

[[ -z $1  ]]&&echo&&head -n 5 $0 | tail -n 1&&echo&&exit
rAw=$1
bNam=`echo -e "${rAw}" | rev | cut -f 1 -d "/" | rev | sed -e "s/[.]fa//g"`

##### List blastn collections #####
printf "Listing subprocess (`date`)"
#ls ../data/00_${bNam}_*_db.fa > ${bNam}_iDx.csv
printf " -- Done (`date`)\n"

##### Form child script & execute #####
printf "Assembling run-script 'hpc_c.sh' (`date`)"
cat hpcHead.sh > ${bNam}_hpc_c.sh
p0=`echo -e "$ SLURM_ARRAY_TASK_ID" | sed -e "s/ //"`
echo -e "#SBATCH -J ${bNam}\n#SBATCH --array=1-`wc -l < ${bNam}_iDx.csv`\n\nbash ${bNam}_dbSum_c.sh ${p0}" >> ${bNam}_hpc_c.sh
printf " -- Done (`date`)\n"

printf "Modifying child-script 'dbSum_c.sh' to ${bNam}_dbSum_c.sh (`date`)"
sed -e "s/iDx/${bNam}_iDx/g" dbSum_c.sh > ${bNam}_dbSum_c.sh
printf " -- Done (`date`)\n"
sbatch ${bNam}_hpc_c.sh
exit
