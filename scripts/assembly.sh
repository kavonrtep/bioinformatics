#PBS -l select=1:ncpus=2:mem=2gb -l walltime=2:00:00
# the line above species that we request one computer node with
# 2 processors(cpu) and 2 GB of RAM. computation will take max 2 hrs

module add trinity-2.6.5
Trinity --max_memory 2G --seqType fq --left WT_01_R1.fastq  --right WT_01_R2.fastq --output output_directory_trinity --CPU 2
