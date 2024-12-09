#!/bin/bash
# author: ph-u
# script: blastn_c.sh
# desc: child script calling blastn sequence alignment
# in: bash blastn_c.sh [$SLURM_ARRAY_TASK_ID]
# out: res/[gene]-[clinical].txt
# arg: 1
# date: 20231216, 20241119

##### Requirements #####
# iDx.csv # List of gene names equivalent to reference gene fasta filename
# iDx.txt # List of NCBI blastn database(s)
########################

tID=$1
i2=`head -n ${tID} iDx.csv | tail -n 1 | cut -f 1 -d ","` # source FASTA file name
i3=`echo -e "${i2}" | rev | cut -f 1 -d "_" | rev` # gene name
i5=`echo -e "${i2}" | rev | sed -e "s/_/@/" -e "s/[|]/@/g" | rev | cut -f 1 -d "@"` # strain name
dbNum=`wc -l < iDx.txt`
for i4 in `seq 1 ${dbNum}`;do
    i1=`head -n ${i4} iDx.txt | tail -n 1`
    echo -e "`date` - ${i1}: ${i4} / ${dbNum}"
## blast
    ./blastn -query ../data/${i2}.fa -db `echo -e "${i1}" | rev | cut -f 1 -d "/" | rev` -out ../data/${i5}_${i3}.txt -outfmt "6 delim=@ qseqid sseqid sstart send sstrand sseq" -max_target_seqs 1000 -subject_besthit
    Rscript blastn2faByGene.r ../data/${i5}_${i3}.txt iDx.txt
done
exit
