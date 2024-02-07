# dNdS gene-by-gene comparative genomics pipeline

## Computational requirements

1. Install ```blastn``` program from NCBI, program within `dNdSGenome` directory [NCBI blast](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
0. Use Genbank annotated database (GCA_*), can be downloaded via the NCBI ```datasets``` toolbox [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
0. Construct NCBI databases, one genome one database, all files (including the raw genomic file (with extension `.fna`) store within `dNdSGenome` directory (use the ```makeblastdb``` tool within the "NCBI blast" toolbox)
