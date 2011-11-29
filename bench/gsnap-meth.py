"""
runner for gsnap on methylation data.
gmap_build, gsnap and samtools must be on the PATH of the calling environment
"""

"""
gmap_build -k 12 -D /usr/local/src/methylcoder/bench/reference -d thaliana_v10 /usr/local/src/methylcoder/bench/reference/thaliana_v10.fasta > /usr/local/src/methylcoder/bench/methylcoder_gsnap/gmap_build.log && 
cmetindex -d thaliana_v10 -F /usr/local/src/methylcoder/bench/reference > gmap_cmetindex.log 2> gmap_cmetindex.error.log

gsnap --quiet-if-excessive -A sam -k 12  --nofails --nthreads 2 -D /usr/local/src/methylcoder/bench/reference --quiet-if-excessive --npaths 1 -d thaliana_v10 --mode=cmet-nonstranded /usr/local/src/methylcoder/bench/reads/WT_endosperm_BS_seq_raw_batch-2.1.fastq.trim /usr/local/src/methylcoder/bench/reads/WT_endosperm_BS_seq_raw_batch-2.2.fastq.trim > /usr/local/src/methylcoder/bench/methylcoder_gsnap/methylcoded.gsnap.sam 2> /usr/local/src/methylcoder/bench/methylcoder_gsnap/gsnap_run.log

gsnap --quiet-if-excessive -A sam -k 12  --nofails --nthreads 2 -D /usr/local/src/methylcoder/bench/reference --quiet-if-excessive --npaths 1 -d thaliana_v10 --mode=cmet-nonstranded /usr/local/src/methylcoder/bench/reads/WT_endosperm_BS_seq_raw_batch-2.1.fastq.trim
"""

import argparse
import sys
import os
import os.path as op
import subprocess as sp

def sh(cmd):
    return sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

def gsnap_meth(reference, reads, out_dir, kmer=15, stranded=False, extra_args=""):
    ref_dir = op.dirname(reference)
    ref_base = op.splitext(op.basename(reference))[0] # locals
    ref_name = op.basename(reference) # locals.
    assert os.access(reference, os.R_OK), ("reference not found / readable")
    assert os.access(ref_dir, os.W_OK), ("%s must be writable" % (ref_dir))

    cmd = "gmap_build -k %(kmer)i -D %(ref_dir)s -d %(ref_base)s %(reference)s"
    cmd %= locals()
    print >>sys.stderr, cmd + "\n"

    cmd_cmet = "cmetindex -d %(ref_base)s -F %(ref_dir)s -k %(kmer)d\n"
    cmd_cmet %= locals()
    print >>sys.stderr, cmd_cmet

    threads = 2 # locals
    mode = ["cmet-nonstranded", "cmet-stranded"][int(stranded)]
    reads_str = " ".join(reads)
    cmd_gsnap = "gsnap --npaths 1 --quiet-if-excessive --nthreads %(threads)s \
        -A sam -k %(kmer)d -D %(ref_dir)s -d %(ref_base)s --mode %(mode)s \
           %(reads_str)s"
    cmd_gsnap += "| samtools view -bSF 4 - > %(out_dir)s/gsnap-meth.u.bam"
    cmd_gsnap %= locals()
    print >>sys.stderr, cmd_gsnap

    cmd_sort = "samtools sort -m 400000 %(out_dir)s/gsnap-meth.u.bam \
            %(out_dir)s/gsnap-meth\n" % locals()
    print >>sys.stderr, cmd_sort

    cmd_index = "samtools index %(out_dir)s/gsnap-meth.bam\n" % locals()
    print >>sys.stderr, cmd_index

    cmd_pileup = "samtools mpileup -f %(reference)s -ABI \
            %(out_dir)s/gsnap-meth.bam\n" % locals()
    print >>sys.stderr, cmd_pileup

    conversion = {"C": "T", "G": "A"}
    if stranded: conversion.update({ "T": "C", "A": "G"})

    print "#chrom\tpos1\tn_same\tn_converted"
    for toks in (l.rstrip("\r\n").split("\t") for l in sh(cmd_pileup).stdout):
        chrom, pos1, ref, coverage, bases, quals = toks
        if coverage == '0': continue
        bases, ref = bases.upper(), ref.upper()

        converted = conversion.get(ref)
        if converted is None: continue
        # TODO: reason out the stranded stuffs

        # . == same on + strand, , == same on - strand
        n_same = sum(1 for b in bases if b == ",")
        # converted is "ACGT" for + strand "acgt" for - strand
        n_converted = sum(1 for b in bases if b == converted)
        print "\t".join((chrom, pos1, str(n_same), str(n_converted)))

        if n_same + n_converted == 0:
            print toks


def run(args):

    if not op.exists(args.out_dir): os.makedirs(args.out_dir)
    reference, kmer = op.abspath(args.reference), args.kmer
    reads = map(op.abspath, args.reads)
    gsnap_meth(reference, reads, args.out_dir, args.kmer, args.stranded,
               args.extra_args)






def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--reference", help="reference fasta")
    p.add_argument("-k", dest="kmer", help="kmer length for gsnap. default:" \
                   + "%(default)i", type=int, default=15)
    p.add_argument("--stranded", action="store_true", default=False, help=\
                   "by default, non-stranded library prep is assumed")
    p.add_argument("--out-dir", dest="out_dir", help="path to output directory",
                  default="gsnap-meth")
    p.add_argument("--extra-args", dest="extra_args", help="extra arguments"
                " to send to gsnap")

    p.add_argument('reads', nargs='+', help='reads files (if 2 files are'
                   'they are assumed to be paired end')

    try:
        args = p.parse_args()
    except:
        p.print_help()
        raise

    if (not len(args.reads) in (1, 2)) or args.reference is None:
        sys.exit(not p.print_help())

    run(args)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
