# dNdS gene-by-gene comparative genomics pipeline

## Computational requirements

1. Install blastn program from NCBI, program within `dNdSGenome/` directory ([NCBI blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html))
0. Use Genbank annotated database (GCA\_\*), can be downloaded via the NCBI datasets toolbox ([NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/))
0. (<=v2.0) Construct NCBI databases, one genome one database, all files (including the raw genomic file (with extension `.fna`) store within `dNdSGenome/` directory (use the `makeblastdb` tool within the "NCBI blast" toolbox)
0. (>=v3.0) Construct one NCBI database with all draft genomes and retain all raw unannotated genomes within the directory tree as the `binHPC2/` directory

## Steps

### Preprocessing

1. `bash 00\_ffn2fa.sh [../relative/path/2/reference\_cds\_genome.ffn]`: require a `raw/` directory at the same level as the `src/` directory
0. `bash 01\_genomeFlanking.sh [../relative/2/reference\_cds\_genome.fa] [../path/2/reference_genomic.fa] [(optional) flanking region length]`: Both `.fa` fasta files should be pre-processed in the above step, converting multiple-lined sequences to single-lined sequences
0. `bash 00\_mkdbBLASTN.sh [../relative/path/2/sourceFasta] [../relative/path/2/dbDirectory]`
0. (<=v2.0) move all generated files fropm the last step, also their respective source fasta files (should be with extension `.fna`), to the `dNdSGenome/` (v1.0) or `binHPC/` (v2.0) directory

### Pipeline (in SLURM-mediated high performance computer)

v3.0 binHPC2/ (doi:)  
v2.0 binHPC/ (doi:)  
v1.0 dNdSGenome/ (doi:)
