mkdir ngs_assembly2
cd ngs_assembly2
# get paired end reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022852/SRR022852_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022852/SRR022852_2.fastq.gz
zcat SRR022852_1.fastq.gz | fastx_trimmer -f 2 -o SRR022852_1_trimmed.fastq
zcat SRR022852_2.fastq.gz | fastx_trimmer -f 2 -o SRR022852_2_trimmed.fastq
velveth run_25_paired 25 -fastq -shortPaired -separate SRR022852_1_trimmed.fastq SRR022852_2_trimmed.fastq 
velvetg run_25_paired -ins_length 350
