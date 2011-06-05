MethylCoder
===========

MethylCoder is a single program that takes of bisulfite-treated reads and
outputs per-base methylation data. It also includes scripts for analysis
and visualization.
In addition to a binary output and a SAM alignment file, the direct output
of methylcoder is a text file that looks like ::

    #seqid  mt  bp  c   t
    1   3   1354    0   1
    1   3   1358    0   1
    1   3   1393    0   1
    1   3   1394    0   1
    1   3   1402    0   1
    1   3   1409    0   1

where columns are reference `seqid` methylation context (type) basepair
location(bp) number of reads where a (c)ytosine was unconverted, number
of reads where where a cytosine was converted to (t)hymine. Making methylation
at every methylable basepair easily calculated as c / (c + t).

.. contents ::

About
=====

This software is developed in the `Fischer Lab`_ . At UC Berkeley.
Please report any requests, bugs, patches, problems, docs to bpederse@gmail.com

It is distributed under the `New BSD License <http://github.com/brentp/methylcode/blob/master/LICENSE>`_


Requirements
============

Python
------

Python 2.6 must be installed along with the following modules.
all of these are available from pypi and as such are installable via
::

  $  easy_install [module]

* `numpy`_ to handle arrays and binary data in python
* `pyfasta`_ to access/index fasta files

matplotlib is required to plot the per-chromosome methylation levels.

C
-

* `bowtie`_ to align the reads to the genome.
* (optional) `gsnap`_ (>= 2011-03-28) alternative aligner. (part of gmap).
* (optional) `sam-tools`_ to view the alignments and processing the reads

Installation
------------
once the above python and c libraries are installed, download methylcoder from:

    http://github.com/brentp/methylcode/tarball/master (tar ball)

and untar; or clone the repository via git::

    git clone git://github.com/brentp/methylcode.git


Then, from the methylcode directory, it is still necessary to run ::

    $ sudo python setup.py install

to install the package into your path. After that, the executable 'methylcoder'
will be available on your path. Running with no arguments will print help.


Input
=====
The input to the pipeline is:

* a reference fasta file with one entry per chromosome in the genome to which
  the reads are to be mapped.
* a fastq  or fasta reads file. all reads must be the same length and must be
  from Eckers/Zilberman bisulfite process (with only 2 possibilities not 4 from
  Cokus protocol).
  If 2 read files are specified, they are assumed to be pair ends and the aligner is
  called appropriately.

Output
======

* a textfile containing columns:
   1) seqid (chromosome)
   2) methylation type (1 to 6 see below)
   3) basepair position (0-indexed)
   4) reads with C at this position
   5) reads with T at this position

  methylation can be calculated as column 4 / (column 5 + column 4)
  4 and 5 are corrected for strand (G, A respectively for - strand).

* a set of binary files for each chromosome in the fasta file. each file
   contains a value for each basepair in the chromosome--many of which will be
   0 if the position is not a C or G. these files contain no headers and can be
   read in any language by specifying the file-type (listed in [square
   brackets] below. these include:

   + methyltype.bin with values between 1 and 6 as described below (value of
     0 means no methylation is possible at this basepair). [encoded as uint8]
   + cs.bin containing the number of reads with C's at each position (same as
     column 4 above). [encoded as uint32]
   + ts.bin containing the number of reads with T's at each position (same as
     column 5 above). [encoded as uint32]

* Methylation type is a value between 1 and 6:
   1) CG  on + strand
   2) CHG on + strand
   3) CHH on + strand
   4) CG  on - strand
   5) CHG on - strand
   6) CHH on - strand

Pipeline
========
You must have:

    1) input reference fasta file to which to align the reads. here: `thaliana_v9.fasta`
    2) a reads file in fastq or fasta format. here: `reads.fastq`.
       if you have paired end reads, they must be specified in order 1, 2.
    3) a directory containing the bowtie and bowtie-build executables.
       (or the path to the gmap/gsnap install directory the gsnap utilities

An example command to run the pipeline is::

    $ methylcoder --bowtie /usr/local/src/bowtie/ \
                  --extra-args "-m 1"
                  --reference /path/to/thaliana_v9.fasta \
                  /path/to/reads.fastq

or using the gsnap aligner on paired-end reads.::

    $ methylcoder --gsnap /usr/local/bin/ \
                  --reference /path/to/thaliana_v9.fasta \
                  /path/to/reads_1.fastq /path/to/reads_2.fastq

Where you must adjust `/path/to/reads.fastq` to point to your BS-treated reads.
This will create the files specified in `Output`_ above, sending the text to
`path/to/reads_methylcoder/methy-data-DATE.txt` where DATE is the current date.
The binary files will be sent to, that same directory as:
`thaliana_v9.fasta.[CHR].methyl.bin` where [CHR] is substituted by each
chromosome in the fasta file. Once bowtie is run once, its output is not
deleted, and methylcoder.py will only re-run bowtie if its input has been
modified since it was run last. *NOTE* if the `methylcoder` executable is
called without any options, it will print help and available command-line
arguments.

Additional args can be sent directly to the aligner as a string to methylcoder.py's
--extra-args parameter. This would look like. ::

    --extra-args "--solexa-quals -k 1 -m 1 --strata"

and that string will be passed directly to the bowtie invocation when it is
called from methylcoder. Whenever 2 fastq files are sent, they are assumed
to be paired-end reads.

Limitations
===========

  + when using bowtie, the reference size must be less than about 2 Gigabases.
    This limitation can be circumvented by splitting the reference into 2 smaller
    reference sequences. For example with human, splitting into 2 fasta files,
    one with chromosomes 1-9 and the other with chromosomes 10+ works well.
    This limitation does not exist when GSNAP is used as the aligner.

Analysis/Visualization
======================

See: http://github.com/brentp/methylcode/wikis/using-samtools-to-view-alignments

Reading
=======
* Eckers paper.
  http://www.nature.com/nature/journal/v462/n7271/extref/nature08514-s1.pdf

* Bowtie Paper:
  Langmead B, Trapnell C, Pop M, Salzberg SL. 2009. Ultrafast and memory-efficient
  alignment of short DNA sequences to the human genome. Genome Biol 10:R25.

* GSNAP paper:
  Wu TD, Nacu S. 2010 Fast and SNP-tolerant detection of complex variants and splicing in short reads.
  Bioinformatics. 26(7):873-81.

Notes
=====

**warning**
when run with bowtie, methylcoder assumes that the Bisulfite converted reads are created using the Zilberman/Ecker method in which BS conversion occurs *after* conversion to solexa library--giving only 2 possibibilities. This is in contrast to the Jacobsen method which gives 4 possiblities. When run with gsnap, it is assumed that the Jacobsen method was used.

.. _`cython`: http://cython.org
.. _`numpy`: http://numpy.scipy.org
.. _`pyfasta`: http://pypi.python.org/pypi/pyfasta/
.. _`h5py`: http://pypi.python.org/pypi/h5py/
.. _`bowtie`: http://bowtie-bio.sourceforge.net/index.shtml
.. _`sam-tools`: http://samtools.sourceforge.net/
.. _`Fischer Lab`: http://epmb.berkeley.edu/facPage/dispFP.php?I=8
.. _`gsnap`: http://research-pub.gene.com/gmap/
