# dNdS gene-by-gene comparative genomics pipeline

## Computational requirements

1. Install `blastn` program from NCBI, program within `dNdSGenome/` directory ([NCBI blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html))
0. Use Genbank annotated database (GCA_*), can be downloaded via the NCBI `datasets` toolbox ([NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/))
0. Construct NCBI databases, one genome one database, all files (including the raw genomic file (with extension `.fna`) store within `dNdSGenome/` directory (use the `makeblastdb` tool within the "NCBI blast" toolbox)

## Steps

### Preprocessing

1. `bash 00_ffn2fa.sh [../relative/path/2/reference_cds_genome.ffn]`: require a `raw/` directory at the same level as the `src/` directory
0. `bash 01_genomeFlanking.sh [../relative/2/reference_cds_genome.fa] [../path/2/reference_genomic.fa] [(optional) flanking region length]`: Both `.fa` fasta files should be pre-processed in the above step, converting multiple-lined sequences to single-lined sequences
0. `bash 00_mkdbBLASTN.sh [../relative/path/2/sourceFasta] [../relative/path/2/dbDirectory]`
0. move all generated files fropm the last step, also their respective source fasta files (should be with extension `.fna`), to the `dNdSGenome/` directory

### Pipeline (in SLURM-mediated high performance computer)

v3.0 `binHPC2/` (doi:)
v2.0 `binHPC/` (doi:)
v1.0 `dNdSGenome/` (doi:)  
