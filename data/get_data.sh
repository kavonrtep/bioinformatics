#!/bin/sh
# download data from NCBI sequence read archive
# example
# fastq-dump --outdir SRA --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID
mkdir -p SRA
# TODO add check if file exists
fastq-dump --outdir SRA --gzip --skip-technical  --readids \
           --read-filter pass --dumpbase --split-3 --maxSpotId 1000 --clip ERX272591
