================
BS-Seq Benchmark
================

There are a number of BS-Seq software implementations available.
This is an attempt to provide a reproducible, fair benchmark to 
aid researcher in choosing their software and the implementors in
creating or improving it.

If you have software you'd like to add, Please contact bpederse@gmail.com
The Benchmark is on paired-end reads, so all the software included
can map paired end reads.

Implementations
===============

 + `MethylCoder`_ 

    - MethylCoder wraps `bowtie`_ and `GSNAP`_ and provides a single
      invocation for the entire mapping process. It outputs SAM and
      per-base methylation counts.

    - Pedersen et al. (2011) MethylCoder: Software Pipeline for Bisulfite-Treated
      Sequences. (in review)


 + `Bismark`_
   
    - Wraps `bowtie`_. Maps to all four possible strand configurations. Aims
      for flexibility. It outputs a custom format and genome-wide summary.
      It has an extra program to calculate per-base methylation.

    - Kreuger F. and Andrews S. (2011) Bismark: A flexible aligner and
      methylation caller for Bisulfite-Seq applications. Bioinformatics,
      27, 1571-1572.


 + `BRAT`_

    - Creates a hash of the reference genome and stores seeds into a hash-table
      against which similarly hashed reads are checked.

    - Harris, E.Y. et al. (2010) BRAT: Bisulfite-treated reads analysis tool.
      Bioinformatics, 26, 572-573.


 + `BSMAP`_

    - Uses a bit-wise masks and a 2-bit encoded genome.

    - Xi, Y. and Li, W. (2009) BSMAP: whole genome bisulfite sequence MAPping program.
      BMC Bioinformatics, 10, 232.


Benchmark
=========

The benchmark is to map 5-Million paired-end reads to the *Arabidopsis thaliana*
genome. The reads are initially 45bp, but they are quality-trimmed before sending
to any of the aligners.

The commands to:
   
 + download the reads and reference 
 + trim the reads
 + install all the above software

is available in `get.sh`_
as new version become available, the can be updated in `params.sh`_

The commands to run the benchmark on all software are available 
in `run.sh`_

.. note:: in order to have the memory-use information reported,
          a recent version of Ubuntu (and likely other distributions)
          is required.


Results
=======

.. note:: these are preliminary results and will be updated as feedback
          is received.

If the program contains multiple steps, the highest memory usage and sum
of the processor time across the steps is reported in the row with 'total'
in the *process* column. In addition, in the *total* rows, the *program*
name appears in bold.

.. note:: the brat data will be added shortly.

====================== ===================== =================== =================== ===================
               program            process         memory (MB)      time (minutes)        pairs-mapped
====================== ===================== =================== =================== ===================
               bismark            bismark                1033               730.7                    
               bismark       bismark.prep                1864                43.2                    
           **bismark**              total                1864               773.8             2035400
                  brat               brat                7254                13.1                    
                  brat         brat.count                2921                 2.1                    
                  brat          brat.trim                   5                 0.3                    
              **brat**              total                7254                15.5             1872932
             **bsmap**              bsmap                2556               349.6             2408400
**methylcoder-bowtie** methylcoder-bowtie                4976               213.5             2029439
 **methylcoder-gsnap**  methylcoder-gsnap                5460               941.3             2649210
====================== ===================== =================== =================== ===================



Suggestions
===========

The intent is that this serve as a valuable resource for those doing methylation
analysis. Please create an issue or send an email to bpederse@gmail.com for any
improvements/additions/suggestions.


.. _`MethylCoder`: https://github.com/brentp/methylcode/
.. _`GSNAP`: http://share.gene.com/gmap/
.. _`bowtie`: http://bowtie-bio.sourceforge.net/
.. _`Bismark`: http://www.bioinformatics.bbsrc.ac.uk/projects/bismark/
.. _`BRAT`: http://compbio.cs.ucr.edu/brat/
.. _`BSMAP`: http://code.google.com/p/bsmap/
.. _`get.sh`: https://github.com/brentp/methylcode/blob/master/bench/get.sh
.. _`params.sh`: https://github.com/brentp/methylcode/blob/master/bench/params.sh
.. _`run.sh`: https://github.com/brentp/methylcode/blob/master/bench/run.sh
