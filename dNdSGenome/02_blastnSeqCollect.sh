#!/bin/bash
# author: ph-u
# script: 02_blastnSeqCollect.sh
# desc: collect blastn matching sequences
# in: bash 02_blastnSeqCollect.sh [../relative/path/2/overall_genes.fa]
# out: NA
# arg: 1
# date: 20240125

[[ -z $1 ]]&&echo&&head -n 5 $0 | tail -n 1&&echo&&exit
rAw=$1
bNam=`echo -e "${rAw}" | rev | cut -f 1 -d "/" | rev | sed -e "s/[.]fa//g"`

##### Segregate fasta to every gene #####
printf "Segregating genes from ${rAw} (`date`)"
grep -n ">" ${rAw} > ../data/${bNam}_meta.txt
printf " -- Done (`date`)\n"

##### Form index file #####
ls ../data/00_*.csv | grep -v "_db" | rev | cut -f 1 -d "/" | rev | sed -e "s/[.]csv$//" -e "s/^00_//" > ${bNam}_iDx.csv

##### Form child script & execute #####
echo -e "Assembling run-script '${bNam}_hpc_c.sh' - `date`"
cat hpcHead.sh > ${bNam}_hpc_c.sh
p0=`echo -e "$ SLURM_ARRAY_TASK_ID" | sed -e "s/ //"`
echo -e "#SBATCH -J ${bNam}\n#SBATCH --array=1-`wc -l < ${bNam}_iDx.csv`\n\nbash ${bNam}_blastn_c.sh ${p0}" >> ${bNam}_hpc_c.sh

printf "Modifying child-script 'seqCollect_c.sh' to ${bNam}_blastn_c.sh (`date`)"
sed -e "s/iDx/${bNam}_iDx/g" seqCollect_c.sh > ${bNam}_blastn_c.sh
printf " -- Done (`date`)\n"
sbatch ${bNam}_hpc_c.sh
exit
