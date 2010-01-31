READS=data/s_1_sequence.txt
RAW=/home/brentp/ssd/s_1_sequence.raw
# convert fastq into raw reads.
# TODO: write a script to get unique reads and
# keep order of original fastq
#sed -n '2,${p;n;n;n;}' $READS > $RAW

python code/run_bowtie.py \
    --bowtie /opt/src/bowtie/bowtie-0.12.1/ \
    --reads $RAW \
    --outdir /home/brentp/ssd/ \
    --mismatches 2 \
    --fasta data/thaliana.fasta \
        > /home/brentp/ssd/s_1.txt
