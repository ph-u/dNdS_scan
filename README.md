# dNdS gene-by-gene comparative genomics pipeline

## Computational requirements

1. Install ```blastn``` program from NCBI, program within `dNdSGenome` directory [NCBI blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
0. Use Genbank annotated database (GCA_*), can be downloaded via the NCBI ```datasets``` toolbox [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
0. Construct NCBI databases, one genome one database, all files (including the raw genomic file (with extension `.fna`) store within `dNdSGenome` directory (use the ```makeblastdb``` tool within the "NCBI blast" toolbox)

## Steps

### Preprocessing

1. `bash 00\_ffn2fa.sh [../relative/path/2/reference\_cds\_genome.ffn]`: require a `raw/` directory at the same level as the `src/` directory
0. `bash 01\_genomeFlanking.sh [../relative/2/reference\_cds\_genome.fa] [../path/2/reference\_genomic.fa] [(optional) flanking region length]`: Both `.fa` fasta files should be pre-processed in the above step, converting multiple-lined sequences to single-lined sequences
0. `bash 00\_mkdbBLASTN.sh [../relative/path/2/sourceFasta] [../relative/path/2/dbDirectory]`
0. move all generated files fropm the last step, also their respective source fasta files (should be with extension `.fna`), to the `dNdSGenome/` directory

### Pipeline (in SLURM-mediated high performance computer)

1. `bash 00\_blastn.sh [../relative/path/2/reference\_cds\_genome.fa]`
0. (optional) `bash 01\_resumeBLASTN.sh [../relative/path/2/reference\_cds\_genome.fa] [last slurm jobID]`
0. (optional) `bash 02\_blastnSeqCollect.sh [../relative/path/2/reference\_cds\_genome.fa]`
0. `bash 03\_varPotential.sh [reference\_cds\_genome]`: the cds file (`reference\_cds\_genome.fa`) name is used as the prefix of downstream data collection files for batch recognition
