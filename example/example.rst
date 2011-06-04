============================
MethylCoder Example Analysis
============================

This example will go through the entire process of acquiring some
reads, trimming and filtering them, sending to MethylCoder and then
two simple downstream analyses.
The two analyses differ in the columns and exact information use
but the overall workflow is similar. This can be used as a template
for many types of analyses.

We will use **Arabidopsis thaliana** as the reference organism, but
this will extend to most organisms.

.. contents ::


Variables
=========

This document will contain everything needed to re-create this analysis.
We will rely on some variables ::

    BOWTIE_VERSION=0.12.7
    REF=reference/thaliana_v10.fasta
    READSA=reads/WT_endosperm
    READSB=reads/WT_embryo
    DIRA=./embryo/
    DIRB=./endosperm/
    DATE=2011-06-03
    BEDTOOLS=/usr/local/src/bedtools/bin/
    GROUP=filo/bin/groupBy


These can be changed to match your organism or reads.
 

Getting The Reference
=====================

We will use TAIR version 10 as the reference genome::

    mkdir -p reference
    rm -f $REF
    for i in `seq 1 5` C M
    do
        wget -O - ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/TAIR10_chr${i}.fas \
                                        >> $REF
    done
    # fix the chromosome names.
    perl -pi -e "s/^>([^\s]+).*/>\1/;tr/C/c/" $REF

That simply downloads the chromosomal Fasta files and puts them in a single reference.
For later analysis, we also want to download the annotations::

    wget -O reference/arabidopsis_thaliana.gff ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
   
And for simplicity, we will convert that to a `BED`_ file with only the chrom, start, end, and name of
the genes.
This analysis can also be done on exons, UTRs or any feature. We will use an awk script to convert
from GFF to BED.::

    # create a simple BED file from the genes in the GFF.
    # make the chr names match those in the reference sequence.
    # convert to 0-based start for BED format.
    awk 'BEGIN { FS="\t"; OFS="\t" } 
        ($3 == "gene"){ 
            chrom=tolower($1)
            if (chrom == "chrm"){ chrom = "mitochondia"; }
            if (chrom == "chrc"){ chrom = "chloroplast"; }
            split($9,gene,"=");
            print chrom,$4 - 1,$5,gene[length(gene)] 
        
        }' reference/arabidopsis_thaliana.gff  > reference/thaliana.genes.bed

The rows in thaliana.genes.bed look like::

    chr1    3630    5899    AT1G01010
    chr1    5927    8737    AT1G01020
    chr1    11648   13714   AT1G01030
    chr1    23145   31227   AT1G01040
    chr1    28499   28706   AT1G01046
    chr1    31169   33153   AT1G01050

.. note:: for many organisms. Comprehensive `BED`_ annotation files are available from the `UCSC`_ genome browser.


Getting The Reads
=================

First we will get the reads for this example we will
compare wild-type embryo tissue to wild-type endosperm::

    mkdir -p reads/
    rm -f ${READSA}_1.fastq ${READSA}_2.fastq
    rm -f ${READSB}_1.fastq ${READSB}_2.fastq

    for i in 1 2 3 4 5
    do
        wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/WT_endosperm_BS_seq_raw_batch-${i}.1.fastq >> ${READSA}_1.fastq
        wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/WT_endosperm_BS_seq_raw_batch-${i}.2.fastq >> ${READSB}_2.fastq

        wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/Embryo_BS_seq_raw_batch-${i}.1.fastq >> ${READSB}_1.fastq
        wget -O - http://dzlab.pmb.berkeley.edu:8080/work/GEO_submission/raw/Embryo_BS_seq_raw_batch-${i}.2.fastq >> ${READSB}_2.fastq

    done

.. note:: This is a *lot* of data, if may be sufficient to download the
   first batch (replace the for loop with `for i in 1`

Trimming and Filtering The Reads
================================

Here, we will trim reads with a quality score below 20 (1 % chance of miscall)
and then discard any reads that are shorter than 32 bases after the trimming.
The script included for trimming keeps reads paired so that if one end of a read
is discarded, the other pair is also discarded. This maintains the concordance
expected by most aligners.

This will require that the `fastx_toolkit`_ is installed. See the section
`Installing Software` for help on installing.

::

    python ../bench/scripts/fastq_pair_filter.py -t 20 -l 32 \
                    ${READSA}_1.fastq ${READSA}_2.fastq
    #
    python ../bench/scripts/fastq_pair_filter.py -t 20 -l 32 \
                    ${READSB}_1.fastq ${READSB}_2.fastq

The resulting files `*.fastq.trim` will be filtered and trimmed.

Aligning with MethylCoder
=========================

Here we send the **trimmed** reads to MethylCoder for alignment. We will use bowtie
as the aligner, but we could as easily use GSNAP.  ::

     
    methylcoder --bowtie bowtie/bowtie-${BOWTIE_VERSION} \
        --outdir $DIRA --extra-args "-m 1 --chunkmbs 256" \
        --mismatches=2 --reference $REF ${READSA}_1.fastq.trim  2> endosperm.log

    methylcoder --bowtie bowtie/bowtie-${BOWTIE_VERSION} \
        --outdir $DIRB --extra-args "-m 1 --chunkmbs 256" \
        --mismatches=2 --reference $REF ${READSB}_1.fastq.trim  2> embryo.log

Here we are allowing 2 mismatches and using only uniquely mapped reads (`-m 1`)
to score methylation.
The output for each of the above includes a summary.
For embryo::

    seqid        total_cs     total_ts     CG       CG_cs        CG_ts        CHG      CHG_cs       CHG_ts       CHH      CHH_cs       CHH_ts
    chloroplast  24526        19199375     0.001601 5192         3237430      0.001356 3685         2712942      0.001180 15649        13249003
    chr1         1843943      39931193     0.174965 1062378      5009561      0.038189 259857       6544629      0.018053 521708       28377003
    chr2         1896314      25126949     0.235495 922927       2996162      0.082545 354926       3944869      0.032889 618461       18185918
    chr3         1952119      30683595     0.211370 1011168      3772719      0.065309 348210       4983560      0.026321 592741       21927316
    chr4         1507677      23520427     0.217731 796152       2860435      0.064232 260695       3797955      0.026040 450830       16862037
    chr5         1974378      35503589     0.196893 1077222      4393893      0.051880 315832       5771894      0.022428 581324       25337802
    mitochondria 91705        8937455      0.023371 37365        1561403      0.015696 25411        1593524      0.004978 28929        5782528
    genome-wide  9290662      182902583    0.170902 4912404      23831603     0.050735 1568616      29349373     0.021200 2809642      129721607

and endosperm::

    seqid        total_cs     total_ts     CG       CG_cs        CG_ts        CHG      CHG_cs       CHG_ts       CHH      CHH_cs       CHH_ts
    chloroplast  11544        8814864      0.001589 2411         1515367      0.001310 1676         1277255      0.001237 7457         6022242
    chr1         1438908      39154125     0.143123 871990       5220599      0.033809 229432       6556769      0.012177 337486       27376757
    chr2         1432483      24465608     0.190348 736847       3134207      0.071958 304455       3926574      0.021981 391181       17404827
    chr3         1513520      30354534     0.170395 823354       4008684      0.057377 306293       5031985      0.017692 383873       21313865
    chr4         1154105      23190205     0.175575 645888       3032816      0.056165 227661       3825789      0.016889 280556       16331600
    chr5         1516515      34924767     0.159267 875830       4623308      0.045367 275518       5797634      0.014684 365167       24503825
    mitochondria 70497        5180495      0.031891 30069        912789       0.023564 22484        931674       0.005350 17944        3336032
    genome-wide  7137572      166084598    0.150804 3986389      22447770     0.047624 1367519      27347680     0.015106 1783664      116289148

Looking at the `CG` column, we can see that the endosperm generally has lower
methylation as was previously reported in this `Science Paper`_
Below we show a more gene-wise analysis.

Downstream Analysis
===================

We will use some linux tools to transform the simple output data to BED format
so we can utilize some common bioinformatics tools. Namely we will used `bedtools`_
to intersect, merge, and group the data.


Differential Methylation
------------------------

Methylation as Bed
******************

For any base, the methylation can be calculated as C / (C + T) or the
proportion of Cytosines that were not converted for the reads covering that
base. We will choose an arbitrary cutoff of 0.4 and say that Cytosines with
a methylation value above that are methylated and those with a value below
that are not methylated. We will write the result to a bed file with the format
chromosome, start, stop, methylated. Where methylated is 0 or 1. ::

    # write to a BED file. value is 1 if methylation is above cutoff or zero otherwise. 
    CUTOFF=0.4
    grep -v '#' $DIRA/methyl-data-${DATE}.txt | awk -v c=$CUTOFF  'BEGIN { OFS="\t" } \
                  { meth=$4/($4 + $5); print $1,$3,$3+1,(meth > c) ? 1 : 0 }' > $DIRA/methylation.bed

    grep -v '#' $DIRB/methyl-data-${DATE}.txt | awk -v c=$CUTOFF  'BEGIN { OFS="\t" } \
                  { meth=$4/($4 + $5); print $1,$3,$3+1,(meth > c) ? 1 : 0 }' > $DIRB/methylation.bed

Those methylation.bed files will now contain rows like: ::

    chr1    306     307     1
    chr1    309     310     1
    chr1    310     311     1
    chr1    313     314     0
    chr1    316     317     0
    chr1    321     322     0
    chr1    322     323     0

here the final columm indicates wether the base described by the
first 3 columns is methylated according to our cutoff. 


Intersect with Genes
********************

We will now intersect this file with are BED annotation file to
get all methylation data associated with a gene. We will do this 
by chromosome to reduce memory usage since the methylation.bed files
are quite large with 1 row for every C or G in the reference genome.::

    rm -f $DIRA/methylation-by-genes.bed
    rm -f $DIRB/methylation-by-genes.bed
    for chrom in `seq 1 5`
    do
        grep chr$chrom $DIRA/methylation.bed | \
            $BEDTOOLS/intersectBed -wao -a reference/thaliana.genes.bed \
                                -b stdin | cut -f 1-4,8 >> $DIRA/methylation.by-genes.bed
        grep chr$chrom $DIRB/methylation.bed | \
            $BEDTOOLS/intersectBed -wao -a reference/thaliana.genes.bed \
                                -b stdin | cut -f 1-4,8 >> $DIRB/methylation.by-genes.bed
    done

The resulting `methylation.by-genes.bed` files now contain 1 row for each C or G
that falls within a gene. 

Count by Gene
*************

We want to get the sum of those within each gene so we use mergeBed 
and sum by gene name. Again we go per chromosome for memory considerations::

    for chrom in `seq 1 5`
    do
        grep chr$chrom $DIRA/methylation.by-genes.bed | \
            $BEDTOOLS/mergeBed -d -1 -nms -scores sum -i stdin \
            | awk 'BEGIN { OFS="\t" } { split($4,names,";"); print $1,$2,$3,names[1],int($5) }' \
            | sort -k4,4 >> $DIRA/counts.bed

        grep chr$chrom $DIRB/methylation.by-genes.bed | \
            $BEDTOOLS/mergeBed -d -1 -nms -scores sum -i stdin \
            | awk 'BEGIN { OFS="\t" } { split($4,names,";"); print $1,$2,$3,names[1],int($5) }' \
            | sort -k4,4 >> $DIRB/counts.bed
    done


Now `counts.bed` will look something like::

    chr1    3630    5899    AT1G01010       6
    chr1    5927    8737    AT1G01020       18
    chr1    11648   13714   AT1G01030       2
    chr1    23145   33153   AT1G01040       273
    chr1    33378   37871   AT1G01060       14
    chr1    38751   40944   AT1G01070       1

indicating that AT1G01010 has 6 methylated Cytosines. 
This BED file can be entered in a genome-browser to get a visual idea of the differences.


Join the Groups
***************

Now we compare the counts in the embryo to the counts in the endosperm.
To do so, we create a single file with columns for gene, embryo, and endosperm. 
We use the linux tool, *join* to join the 2 files on the gene name and then cut
out the columns we need.::

    # join on the name, and only grab the count of methylated C's from each file.
    echo "gene      embryo       endosperm" > embryo.endosperm.counts.txt
    join -t"        " -j 4 $DIRA/counts.bed $DIRB/counts.bed | cut -f 1,5,9 >> embryo.endosperm.counts.txt

DESeq
*****

This file is now in a format we can use with the `R`_ package `DESeq`_ which generally used to
find differential expression for RNA-Seq data, but can be applied to any count data.

.. note:: DESeq is best used with biological replicates

We will use this `R`_ script::

    tbl = read.delim('embryo.endosperm.counts.txt', header=TRUE, row.names=1, stringsAsFactors=TRUE)
    library(DESeq)

    cds = newCountDataSet(tbl,  c("embryo", "endosperm"))
    cds = estimateSizeFactors(cds)
    # NOTE, we dont have replication in this case.
    cds = estimateVarianceFunctions(cds, method="blind")

    res = nbinomTest(cds, "embryo", "endosperm")
    # p-adjusted < 0.05
    resvalid = res[!is.na(res$padj),]
    resSig = resvalid[ resvalid$padj < .05, ]
    write.table(resSig, "embryo.endosperm.sig.genes.txt", row.names=F, sep="\t",
                        quote=F)

With biological replicates, it is possible to do more thorough and sophisticated
analyses.

GO Analysis
***********

The resulting table in `embryo.endosperm.sig.genes.txt` contains 191 genes
with an adjusted p-value of less than 0.05. These can be sent to a web service
like `amigo`_ to look for gene ontology enrichment. 
For this set, we find that only `GO:0003825 alpha,alpha-trehalose-phosphate synthase (UDP-forming) activity`
is over-represented and not at a significant level. 
From here, we may want to change the 0.4 cutoff used above, we may try using only exons, or,
we may look in the upstream or 5' UTR region of the gene. For any of those analyses, we will
create a BED file and follow the steps above.
Note that in any of those cases, we can change a single parameter, either CUTOFF
or the reference BED file and re-run the analysis verbatim.

Fisher's Exact
--------------

It is common in BS-Seq to use the fisher's exact test between 2 groups where
the counts of C's and T's are used to form the contingency table.  Here, we will
use the total counts of T's and C's within the gene on the embryo and endosperm
groups used above. 

Counts as Bed
*************

This time, we create our BED file with one column for the C
counts and one for the T counts.::

    grep -v '#' $DIRA/methyl-data-${DATE}.txt | awk -v c=$CUTOFF  'BEGIN { OFS="\t" } { print $1,$3,$3+1,$4,$5 }' > $DIRA/ct.bed
    grep -v '#' $DIRB/methyl-data-${DATE}.txt | awk -v c=$CUTOFF  'BEGIN { OFS="\t" } { print $1,$3,$3+1,$4,$5 }' > $DIRB/ct.bed
 
Now, ct.bed has rows like::

    chr1    33      34      3       2
    chr1    44      45      1       0
    chr1    45      46      1       0
    chr1    46      47      1       0
    chr1    52      53      2       0

With the final 2 columns describing the number of unconverted and converted reads at the base
location indicated by the first 3 columns.

Intersect with Genes
********************

As before, we intersect these counts with the genes. Again, we will use
`bedtools`_ by chromosome to minimize memory use.::


    #intersect for every cytosine that falls within a gene.
    rm -f $DIRA/ct-by-genes.bed
    rm -f $DIRB/ct-by-genes.bed

    for chrom in `seq 1 5`
    do
        grep chr$chrom reference/thaliana.genes.bed > reference/t.bed
        grep chr$chrom $DIRA/ct.bed | \
            $BEDTOOLS/intersectBed -wao -a reference/t.bed \
                                -b stdin | cut -f 1-4,8,9 | grep -v "\-1" >> $DIRA/ct-by-genes.bed
        grep chr$chrom $DIRB/ct.bed | \
            $BEDTOOLS/intersectBed -wao -a reference/t.bed \
                                -b stdin | cut -f 1-4,8,9 | grep -v "\-1" >> $DIRB/ct-by-genes.bed
    done


Count By Gene
*************

Again, we group by the gene. This time we will use the excellent `groupby`_ tool
which allows us to group by the gene and sum both the C and T columns in a single
command::

    $GROUP -i $DIRA/ct-by-genes.bed -g 4 -c 5,6 -o sum,sum | sort -k 1,1 > $DIRA/ct-grouped.bed
    $GROUP -i $DIRB/ct-by-genes.bed -g 4 -c 5,6 -o sum,sum | sort -k 1,1 > $DIRB/ct-grouped.bed

The resulting file will contain rows like::

    AT1G01040       1514    14719
    AT1G01046       19      477
    AT1G01050       52      3295
    AT1G01060       110     8735
    AT1G01070       3       2491

The columns are gene, unconverted counts, converted counts.

Join the Groups
***************

with columns of gene, number of uncoverted (C's) and converted (T's) reads.
We `join` on the first column to get a single file with both embryo and endosperm.::

    echo "gene      endosperm-c     endosperm-t     embryo-c        embryo-t" > endosperm.embryo.ct.txt
    join -t "       " -j 1 $DIRB/ct-grouped.bed $DIRA/ct-grouped.bed >> endosperm.embryo.ct.txt

That file will contain rows like::

    gene    endosperm-c     endosperm-t     embryo-c        embryo-t
    AT1G01010       49      4770    28      4671
    AT1G01020       73      3932    60      4348
    AT1G01030       33      3731    3       3719
    AT1G01040       1335    14662   1514    14719
    AT1G01046       15      301     19      477

Fisher's Exact Test
*******************

From there, it is very simple to run a fisher's exact test on the values in the column.
We will use the `fisher`_ module for python. With that, a possible script looks like::

    from fisher import pvalue
    import sys

    ct_counts = open(sys.argv[1])
    ngenes = sum(1 for _ in ct_counts)
    ct_counts.seek(0)
    header = ct_counts.readline()

    cutoff = 0.0001 / ngenes
    for gene, ac, at, bc, bt in (line.rstrip().split() for line in ct_counts):
        p = pvalue(*map(int, (ac, at, bc, bt)))
        if p.two_tail > cutoff: continue
        print gene, ac, at, bc, bt, p.two_tail

and we call it like::

    $ python run_fisher.py endosperm.embryo.ct.txt > ct.fisher.sig.txt

Note that we correct for multiple testing by counting the number of genes we
will test and divide the p-value by that number. Even with this very stringent
cutoff, we find 537 *differentially-methylated* genes.

GO Analysis
***********

As before, we can take these genes and send to a GO enrichment tool.
For this geneset, we find the following enrichments::

    GO:0005991 trehalose metabolic process   2.49e-03
    GO:0003825 alpha,alpha-trehalose-phosphate synthase (UDP-forming) activity       5.25e-04
    GO:0000049 tRNA binding  1.25e-02

With the final column indicating the p-value of the enrichment. 
We could continue this analysis with other cutoffs in the fisher test or by
looking at promotor methylation. Again, the workflow would remain the same,
only the annotation BED file would differ.

Installing Software
===================

Much of the software used in this example can be installed following the instructions here 
Note that not all of that is needed for this example.

https://github.com/brentp/methylcode/blob/master/bench/get.sh

Bedtools can be installed as::

    wget http://bedtools.googlecode.com/files/BEDTools.v2.12.0.tar.gz
    tar xzvf BEDTools.v2.12.0.tar.gz
    cd BEDTools-Version-2.12.0/ && make

`groupby`_ is part of Aaron Quinlan's `filo` package and can be installed as::

    git clone https://github.com/arq5x/filo.git
    cd filo && make

The `fisher`_ python module can be installed as::

    wget http://pypi.python.org/packages/source/f/fisher/fisher-0.1.4.tar.gz
    tar xzvf fisher-0.1.4.tar.gz && cd fisher-0.1.4 && sudo python setup.py install

The `run_fisher.py` and `de.R` scripts are included in their entirety in this
document, but can also be found in the scripts/ directory.
The script to trim and filter paired-end reads is also included with MethylCoder
in the benchmarks/scripts directory.


.. _`fastx_toolkit`: http://hannonlab.cshl.edu/fastx_toolkit/
.. _`Science Paper`: http://www.sciencemag.org/content/324/5933/1451.full
.. _`bedtools`: https://github.com/arq5x/bedtools/
.. _`R`: http://www.r-project.org/
.. _`DESeq`: http://genomebiology.com/2010/11/10/R106
.. _`UCSC`: http://genome.ucsc.edu/
.. _`amigo`: http://amigo.geneontology.org/cgi-bin/amigo/term_enrichment
.. _`groupby`: https://github.com/arq5x/filo
.. _`fisher`: http://pypi.python.org/pypi/fisher
.. _`BED`: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
