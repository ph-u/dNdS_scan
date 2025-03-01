#!/bin/env bash
# author: ph-u
# script: dnds_runHead.sh
# desc: Run dnds pipeline (template)
# in: bash dNdS_[acc].sh [slurm array index]
# out: NA
# arg: 1
# date: 20250301

#SBATCH -A [Group Account name]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --requeue
#SBATCH -p icelake-himem
#SBATCH --time=12:00:00
