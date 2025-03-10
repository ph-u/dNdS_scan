#!/bin/env bash
# author: ph-u
# script: masterTemplate.sh
# desc: template pipeline, need modifications/segregating sections
# in: bash masterTemplate.sh [stage] [Group-Account-name] [../relative/path/2/refGenomes_list].txt [../relative/path/2/total_blastdb_accession_list].txt [(stage 1 only) Chunk Size]
# out: NA
# arg: 4
# date: 20250224

[[ -z $4 ]] && head -n 5 $0 | tail -n 1 && exit
[[ $1 -gt 3 ]] && sTage=3 || sTage=$1
[[ $1 -lt 1 ]] && sTage=1 || sTage=$1
gpNam=$2; refG=$3; totDB=$4

if [[ ${sTage} -eq 1 ]];then
  [[ -z $5 ]] && head -n 5 $0 | tail -n 1 && exit
  ChunkSize=$5
##### Stage 0: Prepare accession lists #####
# Note that NCBI accession number format: GCF_000... or GCA_000...
# 1. txt file containing accession numbers of reference genomes, one line per number
# 2. csv/txt file containing accession numbers (clinical/environmental isolates) that will become part of the blast database (can be downloaded using NCBI tool `datasets`)
  mkdir -p ../data #&& apptainer pull docker://ghcr.io/ph-u/dnds_scan:latest
#  printf "Assembling configuration toml file (`date`)"
#  echo -e "[memory]\n    limit = 7088373760" > cgroups.toml
#  printf " -- Done (`date`)\n"

##### Stage 1: Reference genome & blastn database segregation #####
  apptainer run --bind ../data:/data dnds_scan_latest.sif ref ${refG} &
  apptainer run --bind ../data:/data dnds_scan_latest.sif dbChunks ${totDB} ${ChunkSize}

##### Stage 2: blastn database preparation #####
elif [[ ${sTage} -eq 2 ]];then
  [[ -f blastdb-download.sh ]]&&rm blastdb-download.sh
  p0=`echo -e "$ {SLURM_ARRAY_TASK_ID}" | sed -e "s/ //"`

  printf "Assembling blastdb-download run-script (`date`)"
  echo -e "#!/bin/env bash\n# author: ph-u (docker container)\n# in: sbatch blastdb-download.sh\n# date: `date`\n#SBATCH -A ${gpNam}\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --mail-type=NONE\n#SBATCH --requeue\n#SBATCH -p icelake-himem\n#SBATCH --time=12:00:00\n#SBATCH -J db_IN\n#SBATCH --array=1-`ls ../data/*_${totDB} | wc -l`\n\napptainer run --bind ${PWD}/../data:/data dnds_scan_latest.sif db ../data/${p0}_${totDB}" >> blastdb-download.sh
  printf " -- Done (`date`)\n"
  sbatch blastdb-download.sh

else
##### Stage 3: Write & Run dN/dS calculations #####
  printf "Assembling dN/dS run-script (`date`)"

  [[ `ls | grep -e "dNdS_G" | grep -e ".sh" | wc -l` -gt 0 ]]&&rm dNdS_G*.sh
  p0=`echo -e "$ {SLURM_ARRAY_TASK_ID}" | sed -e "s/ //"`

  while read -r L;do
    i=`echo -e "${L}" | rev | sed -e "s/ /@/" | cut -f 1 -d "@" | rev` # accession num
    iSt=`grep -n ${i} iDx.csv | head -n 1 | cut -f 1 -d ":"` # first cds in genome
    iEd=`grep -n ${i} iDx.csv | tail -n 1 | cut -f 1 -d ":"` # last cds in genome
    echo -e "#!/bin/env bash\n# author: ph-u (docker container)\n# in: sbatch dNdS_${i}.sh\n# date: `date`\n#SBATCH -A ${gpNam}\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --mail-type=NONE\n#SBATCH --requeue\n#SBATCH -p icelake-himem\n#SBATCH --time=12:00:00\n#SBATCH -J ${i}\n#SBATCH --array=${iSt}-${iEd}\n\napptainer run --bind ${PWD}/../data:/data dnds_scan_latest.sif dnds ${p0}" >> dNdS_${i}.sh # --apply-cgroups cgroups.toml --memory 6760M
  done < ../data/freqSLURM.txt
  [[ -f dNdS_.sh ]]&&rm dNdS_.sh

  printf " -- Done (`date`)\n"
  for i in `ls dNdS_*.sh`;do sbatch ${i};done
fi
exit
