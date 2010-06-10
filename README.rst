MethylCoder
===========

Python code and shell scripts for fast, simple processing of BiSulfite reads
into methylation data. Also includes scripts for analysis and visualization.
In addition to a binary output, the direct output of methylcoder is a text file
that looks like ::

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

This software is distributed under the BSD License and possible because of
the `Fischer Lab`_ . Please report any bugs, patches, problems, docs to
bpederse@gmail.com

It is distributed under the `New BSD License <http://github.com/brentp/methylcode/blob/master/LICENSE>`_


Requirements
============

Python
------

all of these are available from pypi and as such are installable via
::

  $  easy_install [module]

* `numpy`_ to handle arrays and binary data in python
* `pyfasta`_ to access/index fasta files

C
-

* `bowtie`_ to align the reads to the genome.
* (optional) `sam-tools`_ to view the alignments and processing the reads

Installation
------------
once the above python and c libraries are installed, it is still necessary to
run ::

    $ sudo python setup.py install

to install the package into your path. After that, the executable 'methylcoder'
will be available on your path. Running with no arguments will print help.


Input
=====
The input to the pipeline is:

* a reference fasta file with one entry per chromosome in the genome to which
  the reads are to be mapped.
* a fastq reads file. all reads must be the same length and must be from
  Eckers/Zilberman bisulfite process

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
    2) a reads file in fastq format. here: `reads.fastq`.
    3) bowtie built in a directory. here: `/usr/local/src/bowtie/`
    4) an `out/` directory to send the results.

An example command to run the pipeline is::

    $ methylcoder --bowtie=/usr/local/src/bowtie/ \
                  --reads /path/to/reads.fastq \
                  --outdir out/   \
                  --reference /path/to/thaliana_v9.fasta

Where you must adjust `/path/to/` to the appropriate paths and `outdir` must exist.
This will create the files specified in `Output`_ above, sending the text to
`out/methy-data-DATE.txt` where DATE is the current date. The binary files will
be sent to, for example: `out/thaliana_v9.fasta.[CHR].methyl.bin` where [CHR] is
substituted by each chromosome in the fasta file. Once bowtie is run once,
its output is not deleted, and methylcoder.py will only re-run bowtie if its
input has been modified since it was run last. *NOTE* if the `methylcoder`
executable is called without any options, it will print help and available
command-line arguments.
Additional args can be sent directly to bowtie as a string to methylcoder.py's
--bowtie_args parameter. This would look like. ::

    --bowtie_args "--solexa-quals -k 1 -m 1 --strata"

and that string will be passed directly to the bowtie invocation when it is
called from methylcoder.


Analysis/Visualization
======================

See: http://github.com/brentp/methylcode/wikis/using-samtools-to-view-alignments

Reading
=======
* Eckers paper.
  http://www.nature.com/nature/journal/v462/n7271/extref/nature08514-s1.pdf

* Bowtie Paper:
  Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient
  alignment of short DNA sequences to the human genome. Genome Biol 10:R25.

Notes
=====

**warning**
methylcoder assumes that the Bisulfite converted reads are created
using the Zilberman/Ecker method in which BS conversion occurs *after*
conversion to solexa library--giving only 2 possibibilities. This is in
contrast to the Jacobsen method which gives 4 possiblities. (The code in
methylcoder.py could be made to handle the 2 additional possiblities but
does not do so currently)

.. _`cython`: http://cython.org
.. _`numpy`: http://numpy.scipy.org
.. _`pyfasta`: http://pypi.python.org/pypi/pyfasta/
.. _`h5py`: http://pypi.python.org/pypi/h5py/
.. _`bowtie`: http://bowtie-bio.sourceforge.net/index.shtml
.. _`sam-tools`: http://samtools.sourceforge.net/
.. _`Fischer Lab`: http://epmb.berkeley.edu/facPage/dispFP.php?I=8
