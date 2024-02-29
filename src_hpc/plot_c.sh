#!/bin/bash
# author: ph-u
# script: plot_c.sh
# desc: child script plotting metadata with dN/dS analysis result
# in: bash plot_c.sh [$SLURM_ARRAY_TASK_ID]
# out: NA
# arg: 1
# date: 20240229

tID=$1
i2=`head -n ${tID} iDx.csv | tail -n 1 | cut -f 1 -d ","` # source FASTA file name
i3=`echo -e "${i2}" | rev | cut -f 1 -d "_" | rev` # gene name
i5=`echo -e "${i2}" | rev | sed -e "s/_/@/" -e "s/[|]/@/g" | rev | cut -f 1 -d "@"` # strain name

Rscript varTypeMap.r ${i2} ${i5} ${i3}

exit
