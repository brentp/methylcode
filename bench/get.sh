BRAT_VERSION=1.2.2
BISMARK_VERSION=0.5.0
BSMAP_VERSION=2.0

BOWTIE_VERSION=0.12.7
GSNAP_VERSION=2011-03-28

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
tar xzvf gmap-gsnap-${GSNAP_VERSION}.tar.gz
cd gmap-*${GSNAP_VERSION}*
./configure && make -j3
cd ../../


mkdir -p bsmap && cd bsmap
#wget http://bsmap.googlecode.com/files/bsmap-${BSMAP_VERSION}.tgz
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

cd ../../
