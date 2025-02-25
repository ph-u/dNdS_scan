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
singularity pull ghcr.io/ph-u/dnds:latest

##### Stage 1: Reference genome & blastn database preparation #####

#singularity run dnds ref relative_path_2_refGenome.txt

##### Stage 2: Run dN/dS calculations #####

exit
