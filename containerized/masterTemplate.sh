#!/bin/env bash
# author: ph-u
# script: masterTemplate.sh
# desc: template pipeline, need modifications/segregating sections
# in: bash masterTemplate.sh
# out: NA
# arg: 
# date: 20250224

##### Stage 0: Prepare accession lists #####
# Note that NCBI accession number format: GCF_000... or GCA_000...
# 1. txt file containing accession numbers of reference genomes, one line per number
# 2. csv/txt file containing accession numbers (clinical/environmental isolates) that will become part of the blast database (can be downloaded using NCBI tool `datasets`)
mkdir ../data && apptainer pull docker://ghcr.io/ph-u/dnds_scan:latest

##### Stage 1: Reference genome & blastn database preparation #####
apptainer run --bind ${PWD}:/data dnds_scan_latest.sif ref ../src_log/ecoli_refGenome.txt & # [../relative/path/2/refGenomes_list].txt

apptainer run --bind ${PWD}:/data dnds_scan_latest.sif dbChunks ../src_log/ecoli_blastdbACC.txt # [../relative/path/2/total_blastdb_accession_list].txt
for i in `ls ../data/*_ecoli_blastdbACC.txt`;do
    apptainer run --bind ${PWD}:/data dnds_scan_latest.sif db ../data/5_ecoli_blastdbACC.txt & #[../relative/path/2/blastdb_accession_list].txt
done

##### Stage 2: Run dN/dS calculations #####
apptainer run --bind ${PWD}:/data dnds_scan_latest.sif dnds 1 #[line number of ORF in iDx.csv file]

exit
