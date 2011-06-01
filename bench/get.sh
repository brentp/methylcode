#!/bin/bash

source ./params.sh

mkdir -p brat && \
cd brat && \
wget http://compbio.cs.ucr.edu/brat/downloads/brat-${BRAT_VERSION}.tar.gz
tar xzvf brat-${BRAT_VERSION}.tar.gz && \
cd brat-${BRAT_VERSION} && \
make
cd ../../


mkdir -p bsseeker && \
cd bsseeker
wget http://pellegrini.mcdb.ucla.edu/BS_Seeker/BS_Seeker_files/BS_SEEKER.tgz
tar xzvf BS_SEEKER.tgz
cd ../

mkdir -p bismark
cd bismark
wget http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/bismark_v${BISMARK_VERSION}.tar.gz
tar xzvf bismark_v${BISMARK_VERSION}.tar.gz
cd ../

mkdir -p bowtie
cd bowtie
wget -O bowtie-${BOWTIE_VERSION}.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie/${BOWTIE_VERSION}/bowtie-${BOWTIE_VERSION}-src.zip/download
unzip bowtie-${BOWTIE_VERSION}.zip
cd bowtie-${BOWTIE_VERSION}
make -j3
cd ../../

mkdir -p gsnap && cd gsnap
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-${GSNAP_VERSION}.tar.gz
tar xzvf gmap-gsnap-${GSNAP_VERSION}*.tar.gz
cd gmap-*${GSNAP_VERSION}*
./configure --prefix=`pwd` && make -j3 && make install
cd ../../

# fastx for read trimming.
mkdir -p fastx && cd fastx
wget http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.tar.bz2
bunzip2 libgtextutils-0.6.tar.bz2 && tar xvf libgtextutils-0.6.tar
cd libgtextutils-0.6 && ./configure && make && sudo make install && sudo ldconfig
cd ../
sudo apt-get install pkg-config
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit-${FASTX_VERSION}.tar.bz2
bunzip2 fastx_toolkit-${FASTX_VERSION}.tar.bz2
tar xvf fastx_toolkit-${FASTX_VERSION}.tar
# fastx-toolkit doesn't build with gcc-4.5
cd fastx_toolkit-${FASTX_VERSION} && CC=gcc-4.4 ./configure && make && sudo make install
cd ../../



mkdir -p bsmap && cd bsmap
wget http://bsmap.googlecode.com/files/bsmap-${BSMAP_VERSION}.tgz
tar xzvf bsmap-${BSMAP_VERSION}.tgz
cd bsmap-${BSMAP_VERSION} && make -j3
cd ../../

mkdir -p reference && cd reference
rm -f thaliana_v10.fasta
for i in `seq 1 5` C M
do
    wget -O - ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr${i}.fas >> thaliana_v10.fasta
done
perl -pi -e "s/^>([^\s]+).*/>\1/;tr/C/c/" thaliana_v10.fasta

cd ../

mkdir reads/
wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/WT_endosperm_BS_seq_raw_batch-2.1.fastq | head -n 5000000 > reads/WT_endosperm_BS_seq_raw_batch-2.1.fastq
wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/WT_endosperm_BS_seq_raw_batch-2.2.fastq | head -n 5000000 > reads/WT_endosperm_BS_seq_raw_batch-2.2.fastq

python scripts/fastq_pair_filter.py -t 20 -l 32 reads/WT_endosperm_BS_seq_raw_batch-2.1.fastq reads/WT_endosperm_BS_seq_raw_batch-2.2.fastq
