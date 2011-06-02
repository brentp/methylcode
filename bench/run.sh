#!/bin/bash

source ./params.sh

#############
# methylcoder
#############

# methylcoder-bowtie
#-------------------
/usr/bin/time -f "%M %U" python ../methylcoder/__init__.py --bowtie bowtie/bowtie-${BOWTIE_VERSION} \
        --outdir methylcoder_bowtie --extra-args "-m 1 --chunkmbs 256" \
        --mismatches=2 --reference $REF  $R1 $R2 2> methylcoder_bowtie.time

/usr/bin/time -f "%M %U" python ../methylcoder/__init__.py --bowtie bowtie/bowtie-${BOWTIE_VERSION} \
        --outdir methylcoder_bowtie_2 --extra-args "-m 1 --chunkmbs 256" \
        --mismatches=2 --reference $REF  $R1 $R2 2> methylcoder_bowtie_2.time

# methylcoder-gsnap
#------------------
/usr/bin/time -f "%M %U" python ../methylcoder/__init__.py --gsnap gsnap/gmap-${GSNAP_VERSION}*/bin \
        --outdir methylcoder_gsnap --extra-args "--quiet-if-excessive --npaths 1" \
        --mismatches=2 --reference $REF  $R1 $R2 2> methylcoder_gsnap.time

/usr/bin/time -f "%M %U" python ../methylcoder/__init__.py --gsnap gsnap/gmap-${GSNAP_VERSION}*/bin \
        --outdir methylcoder_gsnap_2 --extra-args "--quiet-if-excessive --npaths 1" \
        --mismatches=2 --reference $REF  $R1 $R2 2> methylcoder_gsnap_2.time


##########
# bsseeker
##########

/usr/bin/time -f "%M %U" \
python bsseeker/Preprocessing_genome.py -f $REF -t N \
        -p `pwd`/bowtie/bowtie-${BOWTIE_VERSION} 2> bsseeker.preprocess.time

/usr/bin/time -f "%M %U" \
python bsseeker/BS_Seeker.py -p `pwd`/bowtie/bowtie-${BOWTIE_VERSION} -m 2 \
                   -t N -o bsseeker.output -i $R1 2> bsseeker.map.time

/usr/bin/time -f "%M %U" \
python bsseeker/BSSout2SAM.py -r $REF -f bsseeker.output > bsseeker/output.sam \
 2> bsseeker.sam.time

#######
# bsmap
#######

/usr/bin/time -f "%M %U" \
bsmap/bsmap-${BSMAP_VERSION}/bsmap -a $R1 -b $R1 -d $REF -o bsmap/output.sam \
                        -v 2 -w 1 -p $PROCESSORS -m 0 -x 2000 2> bsmap.time

######
# brat
######

mkdir -p brat_output/

/usr/bin/time -f "%M %U" \
brat/brat-${BRAT_VERSION}/trim -1 $R1 -2 $R2 -P brat_output/reads -L 33 -m 21 2> brat.trim.time
echo "$REF" > brat_output/ref.txt

/usr/bin/time -f "%M %U" \
brat/brat-${BRAT_VERSION}/brat -r brat_output/ref.txt -o brat_output/brat.out \
      -pe -1 brat_output/reads_reads1.txt -2 brat_output/reads_reads2.txt -bs -m 2 2> brat.time

echo "brat_output/brat.out" > brat_output/brat.out.list
/usr/bin/time -f "%M %U" \
brat/brat-${BRAT_VERSION}/remove-dupl -r brat_output/ref.txt -p \
           brat_output/brat.out.list 2> brat.removedup.time

echo "brat_output/brat.out.nodupl" > brat_output/brat.out.list.nodupl
/usr/bin/time -f "%M %U" \
 brat/brat-${BRAT_VERSION}/acgt-count -r brat_output/ref.txt \
          -P brat_output/acgt -p brat_output/brat.out.list.nodupl \
          2> brat.count.time
#########
# bismark
#########

/usr/bin/time -f "%M %U" \
bismark/bismark_v${BISMARK_VERSION}/bismark_genome_preparation --yes \
    --path_to_bowtie `pwd`/bowtie/bowtie-${BOWTIE_VERSION}/ \
    `pwd`/`dirname $REF` 2>bismark_prep.time

/usr/bin/time -f "%M %U" \
bismark/bismark_v${BISMARK_VERSION}/bismark --chunkmbs 256 --fastq -1 $R1 \
     -2 $R2 --path_to_bowtie `pwd`/bowtie/bowtie-${BOWTIE_VERSION}/ \
     `pwd`/`dirname $REF` 2> bismark.time

