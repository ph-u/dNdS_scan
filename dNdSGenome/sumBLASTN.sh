#!/bin/bash
# author: ph-u
# script: sumBLASTN.sh
# desc: summarize blastn results in tabular format
# in: bash sumBLASTN.sh [strain_id] [gene_id] [database_num]
# out: stdout message
# arg: 3
# date: 20231217, 20231221, 20240105

sTrain=$1
gEne=$2
iPcd=`head -n $3 dbList.txt | tail -n 1`
noHit="NoHits%"
#aPcB=36; aPcA=33 # random sequence 25% similarity by chance; what's the optimum similarity?

eQ=""; locusNam="blastnERR%"; dbSt=""; gSt=""; gLen=""; fBef=""; fAft=""; pM=""
if [[ `grep -e "No hits found" ../res/${sTrain}_${gEne}.txt | wc -l` -gt 0 ]];then
    locusNam=${noHit}
elif [[ `grep -e "Score =" ../res/${sTrain}_${gEne}.txt | wc -l` -gt 0 ]];then

## capture blastn result best match
    if [[ `grep -e "Score =" ../res/${sTrain}_${gEne}.txt | wc -l` -lt 2 ]];then
        lN=$(( `grep -n "Lambda" ../res/${sTrain}_${gEne}.txt | head -n 1 | cut -f 1 -d ":"` -1 ))
    else
        lN=$(( `grep -n "Score =" ../res/${sTrain}_${gEne}.txt | head -n 2 | tail -n 1 | cut -f 1 -d ":"` -1 ))
    fi
    lN0=`grep -n "Score =" ../res/${sTrain}_${gEne}.txt | head -n 1 | cut -f 1 -d ":"`
    head -n ${lN} ../res/${sTrain}_${gEne}.txt | tail -n +${lN0} > ../data/00_${sTrain}_${gEne}_t.txt # capture the best match gene content

## Extract matching gene segment from BLASTN result
    grep -e "Sbjct" ../data/00_${sTrain}_${gEne}_t.txt | tr -s " " > ../data/00_${sTrain}_${gEne}_t0.txt
    dbSt=`head -n 1 ../data/00_${sTrain}_${gEne}_t0.txt | tail -n 1 | cut -f 2 -d " "`
    dbEd=`tail -n 1 ../data/00_${sTrain}_${gEne}_t0.txt | cut -f 4 -d " "`
    [[ ${dbSt} -gt ${dbEd} ]]&&pM="minus"||pM="plus"

## Identify whether top match is a gene analogue by flanking regions, accept only if yes
    rfNm=`grep -e "${gEne}" ../data/${sTrain}-flanking.csv` # flanking regions & details
    rfBf=`echo -e "${rfNm}" | cut -f 5 -d ","` # before gene
    rfAt=`echo -e "${rfNm}" | cut -f 6 -d ","` # after gene
    rfBfLen=`echo -e "${rfBf}" | grep -o "[A-Z]" | wc -l` # length before gene
    rfAtLen=`echo -e "${rfAt}" | grep -o "[A-Z]" | wc -l` # length after gene

### extract the multi-line contig containing matched database sequence
    dbGNm=`grep -e ">" ../res/${sTrain}_${gEne}.txt | head -n 1 | cut -f 1 -d " "`
    grep -n ">" ${iPcd} > ../data/00_${sTrain}_${gEne}_t0.txt
    dbGlN=`grep -n "${dbGNm}" ../data/00_${sTrain}_${gEne}_t0.txt | cut -f 1,2 -d ":"` # start line of the multi-line contig
    dbLn0=`echo -e "${dbGlN}" | cut -f 1 -d ":"` # which header line in the list of headers
    dbLn1=`echo -e "${dbGlN}" | cut -f 2 -d ":"` # where is the header line in seq file
    if [[ ${dbLn0} -eq `wc -l < ../data/00_${sTrain}_${gEne}_t0.txt` ]];then
        dbLs0=`wc -l < ${iPcd}`
    else
        dbLs0=$(( `head -n $(( ${dbLn0} +1 )) ../data/00_${sTrain}_${gEne}_t0.txt | tail -n 1 | cut -f 1 -d ":"` -1 )) # last line of the multi-line contig
    fi
    head -n ${dbLs0} ${iPcd} | tail -n +$(( ${dbLn1} +1 )) > ../data/00_${sTrain}_${gEne}_t0.txt
    wcOneLn=$(( `head -n 1 ../data/00_${sTrain}_${gEne}_t0.txt | wc -c` -1 ))
    wcLstLn=$(( `head -n 1 ../data/00_${sTrain}_${gEne}_t0.txt | wc -c` -1 ))
    dbLen=$(( ${wcOneLn}*$(( `wc -l < ../data/00_${sTrain}_${gEne}_t0.txt` -1 )) +${wcLstLn} ))

### extract same flanking region length on clinical isolate, correct sequence direction
    if [[ ${pM} == "plus" ]];then
        fk1=$(( ${dbSt} -1 )); fk2=$(( ${dbEd} +1 ))
    else
        fk1=$(( ${dbEd} -1 )); fk2=$(( ${dbSt} +1 ))
    fi
    [[ ${rfBfLen} -gt ${rfAtLen} ]]&&lEn=${rfBfLen}||lEn=${rfAtLen}
    fk1X=$(( ${fk1}-${lEn} )); fk2X=$(( ${fk2}+${lEn} ))
    [[ ${fk1} -eq 0 ]]&&fk1Q="x"||fk1Q="v" # smaller flank exist?
    [[ ${fk2} -eq ${dbLen} ]]&&fk2Q="x"||fk2Q="v" # larger flank exist?
    [[ ${fk1X} -lt 1 ]]&&fk1X=1 # smaller flank shorter?
    [[ ${fk2X} -gt ${dbLen} ]]&&fk2X=${dbLen} # larger flank shorter?
# --> smaller number flanking region: ${fk1X}-${fk1}
# --> larger number flanking region: ${fk2}-${fk2X}

# extract smaller number flank
    if [[ ${fk1Q} == "x" ]];then
        fk1Sin=""
    else
        head -n $(( `echo -e "${fk1}/${wcOneLn}" | bc` +1 )) ../data/00_${sTrain}_${gEne}_t0.txt | tail -n +`echo -e "${fk1X}/${wcOneLn}" | bc` > ../data/00_${sTrain}_${gEne}_t1.txt
	fk1Mul=`cat ../data/00_${sTrain}_${gEne}_t1.txt`
	cS0=`echo -e "${fk1X}/${wcOneLn}" | bc`
	cSt=`echo -e "${fk1X}-${wcOneLn}*${cS0}" | bc`
	[[ ${cSt} -eq 0 ]]&&cSt=${wcOneLn}
	fk1Sin=`echo "${fk1Mul//[$'\n']}" | cut -c ${cSt}-$(( ${cSt} +${lEn} ))`
    fi
# extract larger number flank
    if [[ ${fk2Q} == "x" ]];then
        fk2Sin=""
    else
        head -n $(( `echo -e "${fk2X}/${wcOneLn}" | bc` +1 )) ../data/00_${sTrain}_${gEne}_t0.txt | tail -n +`echo -e "${fk2}/${wcOneLn}" | bc` > ../data/00_${sTrain}_${gEne}_t1.txt
	fk2Mul=`cat ../data/00_${sTrain}_${gEne}_t1.txt`
	cS0=`echo -e "${fk2}/${wcOneLn}" | bc`
	cSt=`echo -e "${fk2}-${wcOneLn}*${cS0}" | bc`
	[[ ${cSt} -eq 0 ]]&&cSt=${wcOneLn}
	fk2Sin=`echo "${fk2Mul//[$'\n']}" | cut -c ${cSt}-$(( ${cSt} +${lEn} ))`
    fi
    [[ -f ../data/00_${sTrain}_${gEne}_t1.txt ]] && rm ../data/00_${sTrain}_${gEne}_t1.txt
# if minus strand, rev -> nFlipSeq; assign flanks before & after gene
    if [[ ${pM} == "plus" ]];then
        dbBf=${fk1Sin}; dbAt=${fk2Sin}
    else
        if [[ ${fk2Q} != "x" ]];then
            dbBf0=`echo -e ${fk2Sin} | rev`
            dbBf=`Rscript nFlipSeq.r ${dbBf0}`
        fi
        if [[ ${fk1Q} != "x" ]];then
            dbAt0=`echo -e ${fk1Sin} | rev`
            dbAt=`Rscript nFlipSeq.r ${dbAt0}`
        fi
    fi
# shorten ref flanks if dbFlanks are shorter
    if [[ $(( `echo -e "${dbBf}" | wc -c` -1 )) -lt ${lEn} ]];then
        rfBfLen=$(( `echo -e "${rfBf}" | wc -c` -1   ))
	rfBfLen0=$(( `echo -e "${dbBf}" | wc -c` -1  ))
        if [[ ${rfBfLen0} -eq 0 ]];then
            rfBf=""
        else
            rfBf=`echo -e "${rfBf}" | cut -c $(( ${rfBfLen} -${rfBfLen0} +1  ))-${rfBfLen}`
        fi
    fi
    if [[ $(( `echo -e "${dbAt}" | wc -c` -1 )) -lt ${lEn} ]];then
        rfAtLen0=$(( `echo -e "${dbAt}" | wc -c` -1  ))
        if [[ ${rfAtLen0} -eq 0  ]];then
            rfAt=""
        else
            rfAt=`echo -e "${rfAt}" | cut -c 1-${rfAtLen0}`
        fi
    fi
    [[ -f ../data/00_${sTrain}_${gEne}_t1.txt ]] && rm ../data/00_${sTrain}_${gEne}_t1.txt

### Detect flanking sequences similarity from reference flanks
    fBef=`python3 textSimScore.py "${rfBf}" "${dbBf}"`
    fAft=`python3 textSimScore.py "${rfAt}" "${dbAt}"`
#    if [[ ${fBef} -le ${aPcB} ]] || [[ ${fAft} -le ${aPcA} ]];then locusNam=${noHit};fi

## Proceed data extraction if sequence identity confirmed by flanking regions
    if [[ ${locusNam} != ${noHit} ]];then
        gSt=`grep -e "Query" ../data/00_${sTrain}_${gEne}_t.txt | head -n 1 | sed -e "s/\t/ /g" | sed -e "s/  / /g" | cut -f 2 -d " "`
        gLen=`grep -e "Gaps = " ../data/00_${sTrain}_${gEne}_t.txt | cut -f 2 -d "/" | cut -f 1 -d " "`

## check sequence whether identical
        iD=`grep -e "Identities" ../data/00_${sTrain}_${gEne}_t.txt | cut -f 2 -d "=" | cut -f 1 -d "/" | sed -e "s/ //g"`
        [[ $(( ${gLen} - ${iD} )) -eq 0 ]]&&eQ="ID%"
## check indels
        iD=`grep -e "Identities" ../data/00_${sTrain}_${gEne}_t.txt | cut -f 2 -d "," | cut -f 2 -d "=" | cut -f 1 -d "/" | sed -e "s/ //g"`
        [[ ${iD} != 0 ]]&&eQ="INDEL%"

        if [[ `grep -e "[[]locus_tag=" ../res/${sTrain}_${gEne}.txt | wc -l` -gt 0 ]];then
            locusNam=`grep -e "[[]locus_tag=" ../res/${sTrain}_${gEne}.txt | head -n 1 | sed -e "s/locus_tag=/@/" | cut -f 2 -d "@" | sed -e "s/[]]/@/" | cut -f 1 -d "@"`
        else
            locusNam=`grep -e ">" ../res/${sTrain}_${gEne}.txt | head -n 1 | sed -e "s/>//" | cut -f 1 -d "," | cut -f 1 -d " "`
        fi
    fi
    rm ../data/00_${sTrain}_${gEne}_t.txt
    rm ../data/00_${sTrain}_${gEne}_t0.txt
fi

echo -e "${gEne},${iPcd},${eQ}${locusNam},${dbSt},${gSt},${gLen},${fBef},${fAft},${pM}"
exit
