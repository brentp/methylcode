"""
runner for gsnap on methylation data.
gmap_build, gsnap and samtools must be on the PATH of the calling environment
"""

import argparse
import sys
import os
import os.path as op
import subprocess as sp
import collections

def sh(cmd):
    return sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)


def sequence(chrom, pos1, reference, cache=[None]):
    """
    to get context return 2 bases on either side of pos1

    #>>> sequence('chr1', 1, 'reference/thaliana_v10.fasta')
    #'NNCCC'
    #>>> sequence('chr1', 3, 'reference/thaliana_v10.fasta')
    #'CCCTA'
    """
    if cache[0] is None or cache[0][0] != chrom:
        p = sh('samtools faidx %(reference)s %(chrom)s' % locals())
        # first line is header
        seq = "".join(p.stdout.readlines()[1:]).replace("\n", "")
        cache[0] = (chrom, seq)

    chrom, seq = cache[0]
    s = seq[max(pos1 - 3, 0): pos1 + 2].upper()
    if len(s) == 5: return s
    return ("N" * (5 - len(s))) + s

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

    cmd_sort = "samtools sort %(out_dir)s/gsnap-meth.u.bam \
            %(out_dir)s/gsnap-meth\n" % locals()
    print >>sys.stderr, cmd_sort

    cmd_index = "samtools index %(out_dir)s/gsnap-meth.bam\n" % locals()
    print >>sys.stderr, cmd_index

    cmd_pileup = "samtools mpileup -f %(reference)s -ABIQ 0 \
            %(out_dir)s/gsnap-meth.bam\n" % locals()
    print >>sys.stderr, cmd_pileup

    conversion = {"C": "T", "G": "A"}
    print "#chrom\tpos1\tn_same\tn_converted\tcontext"

    summary = {"CG": collections.defaultdict(int),
               "CHG": collections.defaultdict(int),
               "CHH": collections.defaultdict(int)}

    for toks in (l.rstrip("\r\n").split("\t") for l in sh(cmd_pileup).stdout):
        chrom, pos1, ref, coverage, bases, quals = toks
        ref = ref.upper()
        if coverage == '0': continue

        s = sequence(chrom, int(pos1), reference)
        converted = conversion.get(ref)
        if converted is None: continue
        ctx = get_context(s, ref == "C")

        bases = bases.upper()

        # TODO: account for strand here.
        # . == same on + strand, , == same on - strand
        if ref == "C":
            n_same = sum(1 for b in bases if b == ".")
        else:
            continue
            n_same = sum(1 for b in bases if b == ",")
        n_converted = sum(1 for b in bases if b == converted)

        # converted is "ACGT" for + strand "acgt" for - strand

        summary[ctx[:-1]]['same'] += n_same
        summary[ctx[:-1]]['converted'] += n_converted
        print "\t".join((chrom, pos1, str(n_same), str(n_converted), ctx))

    print >>sys.stderr, "CTX\tc\tt\t%"
    for context in summary:
        s = summary[context]
        print >>sys.stderr, "%s\t%i\t%i\t%.6f" % (context, s["same"],
                        s["converted"], s["same"] / float(s["same"] +
                            s["converted"]))

def get_context(seq5, forward):
    """
    >>> get_context('GACTG', True)
    'CHG+'
    """
    if forward:
        assert seq5[2] == "C", seq5
        if seq5[3] == "G": return "CG+"
        if seq5[4] == "G": return "CHG+"
        return "CHH+"
    else: # reverse complement
        assert seq5[2] == "G", seq5
        if seq5[1] == "C": return "CG-"
        if seq5[0] == "C": return "CHG-"
        return "CHH-"

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
