#!/bin/bash
# author: ph-u
# script: 04_varTypeMap.sh
# desc: Variation type Mapping for gene-wise PAO1 dN/dS result
# in: bash 04_varTypeMap.sh [prefix]
# out: NA
# arg: 1
# date: 20240229

[[ -z $1 ]]&& head -n 5 $0 | tail -n 1 && exit
bNam=$1

##### list pipeline subprocess details #####
printf "Listing subprocess (`date`)"
[[ -f ${bNam}_iDx.csv  ]]&&rm ${bNam}_iDx.csv
for i in `ls ../data/00_${bNam}_*_db.csv | grep -v "genom"`;do
    i0=`echo -e ${i} | rev | cut -f 1 -d "/" | rev | sed -e "s/_db//g"`
    echo -e "${i0}," >> ${bNam}_iDx.csv
done
printf " -- Done (`date`)\n"

##### Form child script & execute #####
printf "Assembling run-script 'hpc_c.sh' (`date`)"
cat ../dNdSGenome/hpcHead.sh > ${bNam}_hpc_c.sh
p0=`echo -e "$ SLURM_ARRAY_TASK_ID" | sed -e "s/ //"`
echo -e "#SBATCH -J ${bNam}\n#SBATCH --array=1-`wc -l < ${bNam}_iDx.csv`\n\nbash ${bNam}_plot_c.sh ${p0}" >> ${bNam}_hpc_c.sh
printf " -- Done (`date`)\n"

printf "Modifying child-script 'plot_c.sh' to ${bNam}_plot_c.sh (`date`)"
sed -e "s/iDx/${bNam}_iDx/g" plot_c.sh > ${bNam}_plot_c.sh
printf " -- Done (`date`)\n"
#sbatch ${bNam}_hpc_c.sh
exit
