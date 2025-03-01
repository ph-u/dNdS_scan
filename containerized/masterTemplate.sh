#!/bin/env bash
# author: ph-u
# script: masterTemplate.sh
# desc: template pipeline, need modifications/segregating sections
# in: bash masterTemplate.sh [../relative/path/2/refGenomes_list].txt [../relative/path/2/total_blastdb_accession_list].txt
# out: NA
# arg: 
# date: 20250224

##### Stage 0: Prepare accession lists #####
# Note that NCBI accession number format: GCF_000... or GCA_000...
# 1. txt file containing accession numbers of reference genomes, one line per number
# 2. csv/txt file containing accession numbers (clinical/environmental isolates) that will become part of the blast database (can be downloaded using NCBI tool `datasets`)
mkdir ../data #&& apptainer pull docker://ghcr.io/ph-u/dnds_scan:latest

##### Stage 1: Reference genome & blastn database preparation #####
apptainer run --bind ${PWD}:/data dnds_scan_latest.sif ref $1 & #../src_log/ecoli_refGenome.txt &

apptainer run --bind ${PWD}:/data dnds_scan_latest.sif dbChunks $2 #../src_log/ecoli_blastdbACC.txt
for i in `ls ../data/*_$2`;do
    apptainer run --bind ${PWD}:/data dnds_scan_latest.sif db ../data/${i} & #[../relative/path/2/blastdb_accession_list].txt
done

##### Stage 2: Write & Run dN/dS calculations #####
printf "Assembling run-script (`date`)"

[[ `ls dNdS_G*.sh | wc -l` -gt 0 ]]&&rm dNdS_G*.sh
p0=`echo -e "$ SLURM_ARRAY_TASK_ID" | sed -e "s/ //"`

while read -r L;do
  i=`echo -e "${L}" | rev | sed -e "s/ /@//" | cut -f 1 -d "@" | rev` # accession num
  iSt=`grep -n ${i} iDx.csv | head -n 1 | cut -f 1 -d ":"` # first cds in genome
  iEd=`grep -n ${i} iDx.csv | tail -n 1 | cut -f 1 -d ":"` # last cds in genome
  cat dnds_runHead.sh > dNdS_${i}.sh
  echo -e "#SBATCH -J ${i}\n#SBATCH --array=${iSt}-${iEd}\n\napptainer run --bind ${PWD}:/data dnds_scan_latest.sif dnds ${p0}" >> dNdS_${i}.sh
done < ../data/freqSLURM.txt

printf " -- Done (`date`)\n"

exit
