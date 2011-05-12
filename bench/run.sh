source ./versions.sh

R1=reads/WT_endosperm_BS_seq_raw_batch-2.1.fastq.trim
R2=reads/WT_endosperm_BS_seq_raw_batch-2.2.fastq.trim
REF=reference/thaliana_v10.fasta

#############
# methylcoder
#############

<<FIRST
/usr/bin/time -f "%M %U" methylcoder --bowtie bowtie/bowtie-${BOWTIE_VERSION} \
        --outdir methylcoder_bowtie --extra-args "-m 1" \
        --mismatches=2 --reference $REF  $R1 $R2 2> methylcoder_bowtie.time

/usr/bin/time -f "%M %U" methylcoder --bowtie bowtie/bowtie-${BOWTIE_VERSION} \
        --outdir methylcoder_bowtie_2 --extra-args "-m 1" \
        --mismatches=2 --reference $REF  $R1 $R2 2> methylcoder_bowtie_2.time
FIRST
/usr/bin/time -f "%M %U" methylcoder --gsnap gsnap/gmap-${GSNAP_VERSION}/bin \
        --outdir methylcoder_gsnap --extra-args "--quiet-if-excessive --npaths 1" \
        --mismatches=2 --reference $REF  $R1 $R2 2> methylcoder_gsnap.time

/usr/bin/time -f "%M %U" methylcoder --gsnap gsnap/gmap-${GSNAP_VERSION}/bin \
        --outdir methylcoder_gsnap_2 --extra-args "--quiet-if-excessive --npaths 1" \
        --mismatches=2 --reference $REF  $R1 $R2 2> methylcoder_gsnap_2.time
