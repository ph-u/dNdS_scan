#!/bin/bash
# author: ph-u
# script: 01_resumeBLASTN.sh
# desc: SLURM-optimized blastn alignment pipeline, resume after getting timeout
# in: bash 01_resumeBLASTN.sh [../relative/path/2/overall_genes.fa] [last slurm jobID]
# out: NA
# arg: 1
# date: 20240107

[[ -z $2 ]]&&echo&&head -n 5 $0 | tail -n 1&&echo&&exit
rAw=$1
lstSlurm=$2
bNam=`echo -e "${rAw}" | rev | cut -f 1 -d "/" | rev | sed -e "s/[.]fa//g"`

##### Segregate fasta to every gene #####
printf "Segregating genes from ${rAw} (`date`)"
grep -n ">" ${rAw} > ../data/${bNam}_meta.txt
printf " -- Done (`date`)\n"

##### Form child script & execute #####
printf "Modifying child-script 'blastn_cRE.sh' to ${bNam}_blastn_c.sh (`date`)"
sed -e "s/iDx/${bNam}_iDx/g" -e "s/xxxxxxx/${lstSlurm}/" blastn_cRE.sh > ${bNam}_blastn_c.sh
printf " -- Done (`date`)\n"
sbatch ${bNam}_hpc_c.sh
exit
