# dNdS gene-by-gene comparative genomics pipeline (v1.0)

### Pipeline (in SLURM-mediated high performance computer)

1. `bash 00_blastn.sh [../relative/path/2/reference_cds_genome.fa]`
0. (optional) `bash 01_resumeBLASTN.sh [../relative/path/2/reference_cds_genome.fa] [last slurm jobID]`
0. (optional) `bash 02_blastnSeqCollect.sh [../relative/path/2/reference_cds_genome.fa]`
0. `bash 03_varPotential.sh [reference_cds_genome]`: the cds file (`reference_cds_genome.fa`) name is used as the prefix of downstream data collection files for batch recognition
