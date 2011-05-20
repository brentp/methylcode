source ./versions.sh

#bismark requires reference to end in .fa
/usr/bin/time -f "%M %U" \
echo bismark/bismark_v${BISMARK_VERSION}/bismark_genome_preparation --yes --path_to_bowtie `pwd`/bowtie/bowtie-${BOWTIE_VERSION}/ `pwd`/`dirname $REF`

