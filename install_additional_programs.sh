#/usr/bash
cd
# FASTQC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
unzip fastqc_v0.11.8.zip
chmod 755 ~/FastQC/fastqc
ln -s ~/FastQC/fastqc ~/bin/fastqc
# when asked set symlink to ~/bin/
apt install fastx-toolkit sickle
