MethylCode
==========

Python code and shell scripts for fast, simple processing of BiSulfite reads
into methylation data. Also includes scripts for analysis and visualization.

.. contents ::

Requirements
============

Python
------

all of these are available from pypi and as such are installable via
::

  $  easy_install [module]


* `numpy`_ for handling arrays and binary data in python
* `pyfasta`_ for easy access and slicing of fasta files
* `cython`_ for fast c-extensions for python
* (optional) `h5py`_ for organizing the output

C
-

* `bowtie`_ for aligning the reads to the genome.
* (optional) `sam-tools`_ for view the alignments and processing the reads

Installation
------------
once the above python and c libraries are installed, it is still necessary to
run ::
    
    $ python setup.py build_ext -i

from the code/ directory to build the c-extension to speed up some python code.


Input
=====
The input to the pipeline is:

* a genomic fasta file with one entry per chromosome in the genome to which
  the reads are to be mapped. 
* a reads file. all reads must be the same length and must be from 
  Eckers/Zilberman process

Output
======

* a textfile containing columns:
   1) seqid (chromosome)
   2) methylation type (1 to 6)
   3) basepair position (0-indexed) 
   4) reads with C at this position
   5) reads with T at this position

  methylation can be calculated as 1 - (column 5 / column 4)
  4 and 5 are corrected for strand (can read G, A respectively for - strand).

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
   + methyl.bin containing the proportion of reads which were methylated at
     this basepair == (c / (c + t)). [encoded as float32]


* Methylation type is a value between 1 and 6:
   + 1 == CG  on + strand
   + 2 == CHG on + strand
   + 3 == CHH on + strand
   + 4 == CG  on - strand
   + 5 == CHG on - strand
   + 6 == CHH on - strand

Pipeline
========
given an input fasta file of `arabidopsis_v8.fasta` and a reads file in fastq
format of `reads.fastq`. First, process the reads file into just the sequence 
(currently this pipeline does not regard the quality info in the fastq file
but may be made to do so in the future) with a unix command like::

    $ sed -n '2,${p;n;n;n;}' reads.fastq > reads.raw 

This will create a file containing only the sequence. cd into the directory
containing the methylcoder.py script (from the directory containing this file
cd into the code/ directory) and run the script with the path to the fasta,
the path to the *directory containing the* bowtie executable (which you have already compiled), the reads
file and the directory where the output will be sent.::

    $ cd code/
    $ python methylcoder.py --bowtie=/path/to/bowtie \
                            --reads /path/to/reads.raw \
                            --outdir r/   \
                            thaliana_v8.fasta 

This will create the files specified in `Output`_ above, sending the text to 
`r/methy-data-DATE.txt` where DATE is the current date. The binary files will
be sent to, for example: thaliana_v8.fasta.[CHR].methyl.bin where [CHR] is 
substituted by each chromosome in the fasta file. Once bowtie is run once,
it's output is not deleted, and methylcoder.py will only re-run bowtie if its
input has been modified since it was run last. *NOTE* if the `methylcoder.py`
script is called without any options, it will print help and its available
commandline arguments.

Given that output, one can then do a sanity check on the output by running::

    $ python sanity_check.py -b -f thaliana_v8.fasta r/thaliana_v8.1.methyl.bin

to check the binary file in the directory '/r' was specified when calling
methylcoder.py above. For a text file, the command is::

    $ python sanity_check.py -t -f thaliana_v8.fasta r/methyl-data-DATE.txt

Because that is reading a text file, it will take a couple minutes, but it 
should *never* fail. Once it's certain that the output is sane, one can create
a moving-window average of the methylation data using the moving_window.py
script. For each input .methyl.bin file, it will create 3 output files, 1 for
each methylation type. So, for the 5 arabidopsis chromosomes, to generate the
15 total moving windows for a window-size of 100, run as::

   $ python python moving_window.py -w 100 r/thaliana_v8.*.methyl.bin

the output files for chromosome 5 will look like:
   * r/thaliana_v8.5.CG.w100.bin
   * r/thaliana_v8.5.CHG.w100.bin
   * r/thaliana_v8.5.CHH.w100.bin

these are written as 32 bit floats.


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
methylcoder.py assumes that the Bisulfite converted reads are created
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
