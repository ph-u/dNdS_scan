#!/bin/bash
# author: ph-u
# script: 03_varPotential.sh
# desc: SLURM-optimized amino acid residue variation protential calculation pipeline
# in: bash 03_varPotential.sh [prefix]
# out: res/[prefix]_[locusTag].{csv,pdf}
# arg: 1
# date: 20240119

[[ -z $1 ]] && echo && head -n 5 $0 | tail -n 1 && echo && exit
pF=$1
mkdir -p ../res

##### List pipeline subprocess details #####
echo -e "Listing dN/dS evaluation files - `date`"
ls ../data/01_${pF}_*.csv > ${pF}_iDx.csv

##### Form child script & execute #####
echo -e "Assembling run-script '${pF}_hpc_c.sh' - `date`"
cat hpcHead.sh > ${pF}_hpc_c.sh
p0=`echo -e "$ SLURM_ARRAY_TASK_ID" | sed -e "s/ //"`
echo -e "#SBATCH -J ${pF}\n#SBATCH --array=1-`wc -l < ${pF}_iDx.csv`\n\nbash ${pF}_rDNDS.sh ${p0}" >> ${pF}_hpc_c.sh
echo -e "Modifying schild-script 'rDNDS.sh' to '${pF}_rDNDS.sh' - `date`"
sed -e "s/iDx/${pF}_iDx/g" rDNDS.sh > ${pF}_rDNDS.sh
echo -e "Done - `date`"
sbatch ${pF}_hpc_c.sh
exit
