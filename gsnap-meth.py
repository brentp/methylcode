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

def sh(cmd, log, wait=True):
    print >>sys.stderr, "[running command] %s" % cmd
    p = sp.Popen(cmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    if wait:
        for line in p.stderr:
            print >>log, line,
        p.wait()
        if p.returncode != 0 or "aborted" in p.stderr.read():
            sys.exit(p.returncode)
    return p

def nopen(f, mode="rb"):
    if not isinstance(f, basestring): return f
    return open(f, mode)

def sequence(chrom, pos1, reference, log, cache=[None]):
    """
    to get context return 2 bases on either side of pos1

    #>>> sequence('chr1', 1, 'reference/thaliana_v10.fasta')
    #'NNCCC'
    #>>> sequence('chr1', 3, 'reference/thaliana_v10.fasta')
    #'CCCTA'
    """
    if cache[0] is None or cache[0][0] != chrom:
        p = sh('samtools faidx %(reference)s %(chrom)s' % locals(), log,
                                        wait=False)
        # first line is header
        seq = "".join(p.stdout.readlines()[1:]).replace("\n", "")
        cache[0] = (chrom, seq)

    chrom, seq = cache[0]
    s = seq[max(pos1 - 3, 0): pos1 + 2].upper()
    if len(s) == 5: return s
    return ("N" * (5 - len(s))) + s




def gmap_built(ref_dir, ref_base):
    for ext in ("chromosome", "contig.iit", "ref12153offsetscomp",
                "chromosome.iit", "genomecomp", "ref12153positions",
                "chrsubset", "maps", "version", "contig",
                "ref12153gammaptrs"):

        f = "%s/%s/%s.%s" % (ref_dir, ref_base, ref_base, ext)
        if not op.exists(f):
            return False
    return True

def cmetindexed(ref_dir, ref_base, kmer):

    for ext in ("metct12%i3positions", "metct12%i3offsetscomp",
                "metct12%i3gammaptrs", "metga12%i3offsetscomp",
                "metga12%i3positions", "metga12%i3gammaptrs"):
        ext = ext % kmer
        f = "%s/%s/%s.%s" % (ref_dir, ref_base, ref_base, ext)
        if not op.exists(f):
            return False
    return True

def gsnap_meth(reference, reads, out_dir, kmer=15, stranded=False, extra_args=""):
    ref_dir = op.dirname(reference)
    ref_base = op.splitext(op.basename(reference))[0] # locals
    ref_name = op.basename(reference) # locals.
    assert os.access(reference, os.R_OK), ("reference not found / readable")
    assert os.access(ref_dir, os.W_OK), ("%s must be writable" % (ref_dir))

    threads = 4 # locals
    log = open(out_dir + "/gsnap-meth-run.log", "w")

    print >>sys.stderr, ("writing log to: %s" % log.name)
    cmd = "gmap_build -w 1 -k %(kmer)i -D %(ref_dir)s -d %(ref_base)s %(reference)s"
    cmd %= locals()
    if not gmap_built(ref_dir, ref_base):
        sh(cmd, log)
    else:
        print >>sys.stderr, "[skipping command] %s" % cmd

    cmd_cmet = "cmetindex -d %(ref_base)s -F %(ref_dir)s -k %(kmer)d\n"
    cmd_cmet %= locals()
    if not cmetindexed(ref_dir, ref_base, kmer):
        sh(cmd_cmet, log)
    else:
        print >>sys.stderr, "[skipping command] %s" % cmd

    mode = ["cmet-nonstranded", "cmet-stranded"][int(stranded)]
    reads_str = " ".join(reads)
    cmd_gsnap = "gsnap --npaths 1 --quiet-if-excessive --nthreads %(threads)s \
        -A sam -k %(kmer)d -D %(ref_dir)s -d %(ref_base)s --mode %(mode)s \
         %(extra_args)s  %(reads_str)s"
    cmd_gsnap += "| samtools view -bSF 4 - > %(out_dir)s/gsnap-meth.u.bam"
    cmd_gsnap %= locals()
    sh(cmd_gsnap, log)

    cmd_sort = "samtools sort %(out_dir)s/gsnap-meth.u.bam \
            %(out_dir)s/gsnap-meth\n" % locals()
    sh(cmd_sort, log)

    cmd_index = "samtools index %(out_dir)s/gsnap-meth.bam\n" % locals()
    sh(cmd_index, log)

    summarize_bam("%(out_dir)s/gsnap-meth.bam" % locals(), reference, log)

def summarize_bam(bam, reference, log=sys.stderr):
    assert op.exists(bam)
    assert op.exists(reference)
    cmd_pileup = "samtools mpileup -f %(reference)s -ABIQ 0 \
            %(bam)s\n" % locals()
    summarize_pileup(sh(cmd_pileup, log, wait=False).stdout, reference, log)

def summarize_pileup(fpileup, reference, log):

    conversion = {"C": "T", "G": "a"}
    print "#chrom\tpos1\tn_same\tn_converted\tcontext"

    summary = collections.defaultdict(lambda:
              {"CG": collections.defaultdict(int),
               "CHG": collections.defaultdict(int),
               "CHH": collections.defaultdict(int)})
    for toks in (l.rstrip("\r\n").split("\t") for l in nopen(fpileup)):
        chrom, pos1, ref, coverage, bases, quals = toks
        if coverage == '0': continue
        ref = ref.upper()
        converted = conversion.get(ref)
        if converted is None: continue


        s = sequence(chrom, int(pos1), reference, log)
        ctx = get_context(s, ref == "C")


        # . == same on + strand, , == same on - strand
        if ref == "C":
            n_same = sum(1 for b in bases if b in ".")
        else:
            n_same = sum(1 for b in bases if b in ",")

        n_converted = sum(1 for b in bases if b == converted)

        if n_same + n_converted == 0: continue

        summary[chrom][ctx[:-1]]['same'] += n_same
        summary[chrom][ctx[:-1]]['converted'] += n_converted
        print "\t".join((chrom, pos1, str(n_same), str(n_converted), ctx))

    print >>sys.stderr, "#chrom\tctx\tc\tt\tpct_methylation"
    for chrom in summary:
        for context in summary[chrom]:
            s = summary[chrom][context]
            print >>sys.stderr, "%s\t%s\t%i\t%i\t%.6f" % (chrom, context, s["same"],
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
