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
apptainer pull docker://ghcr.io/ph-u/dnds_scan:latest

##### Stage 1: Reference genome & blastn database preparation #####
apptainer run --bind ${PWD}:/data dnds_scan_latest.sif ref [../relative/path/2/refGenomes_list].txt
apptainer run --bind ${PWD}:/data dnds_scan_latest.sif db [../relative/path/2/blastdb_accession_list].txt

##### Stage 2: Run dN/dS calculations #####
apptainer run --bind ${PWD}:/data dnds_scan_latest.sif dnds [line number of ORF in iDx.csv file] [refAccessionNumber]

exit
