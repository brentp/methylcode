#!/bin/bash

source ./params.sh


##################
# cleanup and prep
##################

# directory to save the timings.
rm -rf bench-results/
mkdir -p bench-results/
rm -f reads/*.c2t* reads/*.txt

mv $REF ./t.fasta
rm -rf reference/*
mv ./t.fasta $REF
rm -rf ./reference_genome # bsseeker.
rm -rf brat_output/
rm -rf methylcoder_bowtie methylcoder_bowtie_2 methylcoder_gsnap methylcoder_gsnap_2

#############
# methylcoder
#############

# methylcoder-bowtie
#-------------------

/usr/bin/time -f "%M %U" python ../methylcoder/__init__.py --bowtie bowtie/bowtie-${BOWTIE_VERSION} \
        --outdir methylcoder_bowtie --extra-args "-m 1 --chunkmbs 256" \
        --mismatches=2 --reference $REF $R1 $R2 2> bench-results/methylcoder-bowtie.time

/usr/bin/time -f "%M %U" python ../methylcoder/__init__.py --bowtie bowtie/bowtie-${BOWTIE_VERSION} \
        --outdir methylcoder_bowtie_2 --extra-args "-m 1 --chunkmbs 256" \
        --mismatches=2 --reference $REF $R1 $R2 2> bench-results/methylcoder-bowtie-existing-index.time

# methylcoder-gsnap
#------------------
/usr/bin/time -f "%M %U" python ../methylcoder/__init__.py --gsnap gsnap/gmap-${GSNAP_VERSION}*/bin \
        --outdir methylcoder_gsnap --extra-args "--quiet-if-excessive --npaths 1" \
        --mismatches=2 --reference $REF $R1 $R2 2> bench-results/methylcoder-gsnap.time

/usr/bin/time -f "%M %U" python ../methylcoder/__init__.py --gsnap gsnap/gmap-${GSNAP_VERSION}*/bin \
        --outdir methylcoder_gsnap_2 --extra-args "--quiet-if-excessive --npaths 1" \
        --mismatches=2 --reference $REF $R1 $R2 2> bench-results/methylcoder-gsnap-existing-index.time


##########
# bsseeker
##########

/usr/bin/time -f "%M %U" \
python bsseeker/Preprocessing_genome.py -f $REF -t N \
        -p `pwd`/bowtie/bowtie-${BOWTIE_VERSION} 2> bench-results/bsseeker.preprocess.time

/usr/bin/time -f "%M %U" \
python bsseeker/BS_Seeker.py -p `pwd`/bowtie/bowtie-${BOWTIE_VERSION} -m 2 \
                   -t N -o bsseeker/bsseeker.output -i $R1 2> bench-results/bsseeker.mapping.time

<<BROKEN
# CANT GET THIS TO WORK
/usr/bin/time -f "%M %U" \
python bsseeker/BSSout2SAM.py -r $REF -f bsseeker/bsseeker.output > bsseeker/output.sam \
 2> bsseeker.to-sam.time
BROKEN

#######
# bsmap
#######

/usr/bin/time -f "%M %U" \
bsmap/bsmap-${BSMAP_VERSION}/bsmap -a $R1 -b $R2 -d $REF -o bsmap/output.sam \
                        -v 2 -w 2 -p $PROCESSORS -m 0 -x 250 2> bench-results/bsmap.time

######
# brat
######

mkdir -p brat_output/

# brat requires each chr in a seperate file.
pyfasta split $REF --header "reference/%(seqid)s.single.fasta"
mkdir -p brat_output/
ls reference/*.single.fasta > brat_output/ref.txt


/usr/bin/time -f "%M %U" \
brat/brat-${BRAT_VERSION}/trim -1 $R1 -2 $R2 -P brat_output/reads -L 33 -m 21 2> bench-results/brat.trim.time

/usr/bin/time -f "%M %U" \
brat/brat-${BRAT_VERSION}/brat -r brat_output/ref.txt -o brat_output/brat.out \
      -i 0 -a 300 -m 2 -pe -bs \
      -1 brat_output/reads_reads1.txt \
      -2 brat_output/reads_reads2.txt 2> bench-results/brat.time

echo "brat_output/brat.out" > brat_output/brat.out.list
/usr/bin/time -f "%M %U" \
 brat/brat-${BRAT_VERSION}/acgt-count -r brat_output/ref.txt \
          -P brat_output/acgt -p brat_output/brat.out.list \
          2> bench-results/brat.count.time


#########
# bismark
#########

/usr/bin/time -f "%M %U" \
bismark/bismark_v${BISMARK_VERSION}/bismark_genome_preparation --yes \
    --path_to_bowtie `pwd`/bowtie/bowtie-${BOWTIE_VERSION}/ \
    `pwd`/`dirname $REF` 2>bench-results/bismark.prep.time

/usr/bin/time -f "%M %U" \
bismark/bismark_v${BISMARK_VERSION}/bismark --chunkmbs 256 --fastq -1 $R1 \
     -2 $R2 --path_to_bowtie `pwd`/bowtie/bowtie-${BOWTIE_VERSION}/ --best \
     `pwd`/`dirname $REF` 2> bench-results/bismark.time

############
#count reads
############


wc -l bsseeker/bsseeker.output | awk '{ print $1 / 2 }' > bench-results/bsseeker.count
wc -l ${R1}_bismark_pe.txt | awk '{ print $0 - 1 }' > bench-results/bismark.count
samtools view -F 0x4 -F 0x100 -S bsmap/output.sam | grep -c "=" | awk '{ print $0 / 2 }' > bench-results/bsmap.count
samtools view -F 0x04 -S methylcoder_gsnap/methylcoded.gsnap.sam | grep -c "=" | awk '{ print $0 / 2 }' > bench-results/methylcoder-gsnap.count
samtools view -c -F 0x4  -S methylcoder_bowtie/methylcoded.sam | awk '{ print $0 / 2 }' > bench-results/methylcoder-bowtie.count

python scripts/make-table-from-times.py bench-results/
