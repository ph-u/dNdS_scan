# Quantifying natural selection

> [!NOTE]
> Wang _et al_ 2025. Nature Communications. DOI: [10.1038/s41467-025-57532-z](https://doi.org/10.1038/s41467-025-57532-z)

## Quick start

1. `[[ -f dnds_scan_latest.sif ]] && rm dnds_scan_latest.sif && apptainer pull docker://ghcr.io/ph-u/dnds_scan:latest && apptainer run dnds_scan_latest.sif`
0. `bash masterTemplate.sh 1 proj-account ref-genome-accession-list.txt NCBI-accession-list.txt 20 &`
0. `bash masterTemplate.sh 2 proj-account ref-genome-accession-list.txt NCBI-accession-list.txt`
0. `bash masterTemplate.sh 3 proj-account ref-genome-accession-list.txt NCBI-accession-list.txt`

> [!WARNING]
> Might need to **run step 2 multiple times** to get all genomes for blast databases due to NCBI refusing frequent connections.

## Program Inputs
- `proj-account`: the SLURM group account you have access to
- `ref-genome-accession-list.txt`: a list of accession numbers (`GCA_*` or `GCF_*`) that are used as reference genome(s)
- `NCBI-accession-list.txt`: a list of accession numbers that will be made into NCBI blast databases; you can use a metadata summary file (with accession numbers as the first column, with/without header) directly downloaded from the NCBI `datasets` tool (CSV format is acceptable)
- `20`: recommended number of genomes in a blast database for whole genome d<sub>N</sub>/d<sub>S</sub> scanning, could be higher for short ORF sequences (highest number tested: 750)

## Computational requirements

1. `apptainer`/`singularity` installed on High Performance Computer cluster
0. SLURM manager - `masterTemplate.sh` need slight modifications for other systems, such as PBS

> [!TIP]
> You may extract the `ref-genome-accession-list.txt` and `NCBI-accession-list.txt` using the [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) tool

## Other Individual commands available

- `apptainer run --bind ../data:/data dnds_scan_latest.sif ref ref-genome-accession-list.txt`: format reference genome only
- `apptainer run --bind ../data:/data dnds_scan_latest.sif db NCBI-accession-list-subset.txt`: download and generate one blastdb using the listed accession numbers from the file provided
- `apptainer run --bind ../data:/data dnds_scan_latest.sif dbChunks NCBI-accession-list.txt 20`: segregate the full list of accession number into subsets of x accessions per group; in this command x=20
- `apptainer run --bind ../data:/data dnds_scan_latest.sif dnds 2`: run one d<sub>N</sub>/d<sub>S</sub> scan on one gene, indicated by its line number (line = 2) in the overall gene list `iDx.csv`
- `apptainer run --bind ../data:/data dnds_scan_latest.sif reCon.r ../data/[xxx]--rDNDS.csv`: run the residue-level d<sub>N</sub>/d<sub>S</sub> reconstruction script by providing a rolling d<sub>N</sub>/d<sub>S</sub> output from the program
- `apptainer shell --bind ../data:/data dnds_scan_latest.sif`: interactive shell
