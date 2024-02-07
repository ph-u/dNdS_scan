#!/bin/env Rscript
# author: ph-u
# script: listSum.r
# desc: summary of scan result
# in: Rscript listSum.r
# out: stdout message
# arg: 0
# dtae: 20240129

for(i in list.files("../data","_db.csv")){
    a0 = read.csv(paste0("../data/",i), header = T)
    print(paste(i, ";", max(a0$seq2), "/ 428 ;", round(max(a0$seq2)/428*100,0),"%"))
    print(round(table(a0$rAnge)/nrow(a0)*100,0))
    print(table(a0$rAnge))
};rm(i,a0)
