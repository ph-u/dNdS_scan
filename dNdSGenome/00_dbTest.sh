#!/bin/bash
# author: ph-u
# script: 00_dbTest.sh
# desc: test blastn database integrity
# in: bash 00_dbTest.sh
# out: [stdout]
# arg: 0
# date: 20231222

for i in `ls *.ndb`;do
    echo ${i}
    blastn -query ../data/g_raw_PA2683.1.fa -db `echo -e ${i} | sed -e "s/[.]ndb//"` -out t.txt
done
rm t.txt
exit
