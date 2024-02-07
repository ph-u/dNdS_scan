#!/bin/bash
# author: ph-u
# script: hpcHead.sh
# desc: get reassemble & call respective script to run
# in: bash hpc_c.sh [slurm array index]
# out: NA
# arg: 1
# date: 20231216

#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --requeue
#SBATCH -p icelake-himem
#SBATCH --time=12:00:00
