#!/bin/nash
sudo apt-get install build-essential libz-dev libncurses5-dev jellyfish bowtie ncbi-blast+ python3-numpy python-numpy cmake

cd ~/Downloads
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2
cd samtools-1.3.1
make

sudo make install
cd ~/Downloads

wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz
tar xzfv Salmon-0.9.1_linux_x86_64.tar.gz
cd Salmon-latest_linux_x86_64/bin
sudo ln -s $PWD/salmon /usr/local/bin/salmon
cd ~/Downloads
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.4.tar.gz
tar zxfv Trinity-v2.8.4.tar.gz
cd trinityrnaseq-Trinity-v2.8.4/
make
sudo ln -s  $PWD/Trinity /usr/local/bin/
