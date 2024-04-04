#!/bin/bash
# author: ph-u
# script: dbSum_c.sh
# desc: schild script for calling sequence summary
# in: bash dbSum_c.sh [$SLURM_ARRAY_TASK_ID]
# out: data/[bNam]_[gene]--dbSum.csv
# arg: 1
# date: 20240331

tID=$1
i0=`head -n ${tID} iDx.csv | tail -n 1 | rev | cut -f 1 -d "/" | rev`
Rscript dbSum.r ${i0}
Rscript rDNDS.r ${i0}
Rscript reCon.r ${i0}
exit
