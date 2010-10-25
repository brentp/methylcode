"""
compare 2 sets of data that have been run through the methylcode pipeline.
usage:

    % prog [options] dir_a/ dir_b/

"""
from fisher import pvalue
import sys
import os.path as op
sys.path.insert(0, op.join(op.dirname(__file__), "../methylcoder"))
from methyl import MethylGroup
from flatfeature import Flat

def bin_setup(chr, adir, bdir, context):
    am = MethylGroup(adir)[chr]
    bm = MethylGroup(bdir)[chr]
    return am.as_context(context), bm.as_context(context)

def run_compare(flat, adir, bdir, context, window, binary, pvalue_cutoff, ratio_range):
    #fh = open('fisher.different.%s.%ibp.gff' % (context, window), 'w')
    fh = sys.stdout
    print >>sys.stderr, "writing to:", fh.name
    print >>fh, "##gff-version 3"
    for chr in flat.seqids:
        try:
            bp_max = len(flat.fasta[chr])
        except KeyError:
            print >>sys.stderr, chr, "not found. skipping"
            continue
        (a_cs, a_ts, a_mask), (b_cs, b_ts, b_mask) = bin_setup(chr, adir, bdir, context)
        for start in xrange(0, bp_max + window, window):
            end = min(start + window, bp_max)
            if start == end: continue
            a_t_count = a_ts[start:end].sum()
            a_c_count = a_cs[start:end].sum()
            b_t_count = b_ts[start:end].sum()
            b_c_count = b_cs[start:end].sum()

            p = pvalue(a_t_count, a_c_count, b_t_count, b_c_count)
            pv = float(p.two_tail)

            if not binary and pv > pvalue_cutoff: continue
            gc = f.fasta[chr][start:end].upper()
            gc = gc.count("G") + gc.count("C")

            # if a_tot or b_tot == 0, then use 'na'
            a_tot = float(a_c_count + a_t_count)
            a_methyl = (a_c_count / a_tot) if a_tot != 0 else 0#None

            b_tot = float(b_c_count + b_t_count)
            b_methyl = (b_c_count / b_tot) if b_tot !=0 else 0#None
            #strand = "+" if a_methyl > b_methyl else "-"
            strand = "."
            # TODO: use absolute?
            plot = a_methyl - b_methyl if not None in (a_methyl, b_methyl) else 'na'

            # scale by total.
            plot = plot / (a_methyl + b_methyl)

            #print plot, a_methyl, b_methyl
            #if plot == 'na': continue
            if binary:
                if plot != 'na':
                    plot == 1 if (ratio_range[0] <= plot <= ratio_range[1]) else 0
            else:
                if not (ratio_range[0] <= plot <= ratio_range[1]):
                    #print >>sys.stderr, "skipping because of ratio range."
                    continue

            if binary and plot != 'na': plot = 0 if pv > pvalue_cutoff else 1

            attrs="p=%.3G;ac=%i;at=%i;bc=%i;bt=%i;gc=%i;plot=%.3G" % \
                        (pv, a_c_count, a_t_count, b_c_count, b_t_count, gc, plot)
            accns = flat.get_features_in_region(chr, start + 1, end)
            accns = [a["accn"] for a in accns]
            if accns:
                attrs +=";accns=" + ",".join(accns)
            print >>fh, "\t".join(map(str, [chr, "methylation", "dmc", start + 1, end, plot, strand, ".", attrs]))

if __name__ == "__main__":
    import optparse, sys
    p = optparse.OptionParser(__doc__)
    p.add_option("--fasta", dest="fasta", help="path to fasta",
                 default="/labdata/thaliana_v9/thaliana_v9.fasta")
    p.add_option("--flat", dest="flat", help="path to flat file",
                 default="/labdata/thaliana_v9/thaliana_v9.flat")
    p.add_option("--window", dest="window", help="window size",
                 type='int', default=50)
    p.add_option("--context", dest="context", help="methylation context "
                 "one of CG, CHG, or CHH")
    p.add_option("-b", dest="binary", action="store_true", default=False,
                 help="use either 1 or 0 for ratio")
    p.add_option("-p", dest="pvalue", type="float", default=1.0,
                 help="p value cutoff")
    p.add_option("--ratio_range", dest="ratio_range", help="optional: range of ratios"
                 " to include. specify in format low:high. e.g.: 0.2:0.8")
    opts, args = p.parse_args()
    if len(args) != 2 or not (opts.context in ('CG', 'CHG', 'CHH')):
        sys.exit(not p.print_help())
    f = Flat(opts.flat, opts.fasta)
    dir_a, dir_b = args
    assert op.exists(dir_a)
    assert op.exists(dir_b)
    rr = (-2, 2)
    if opts.ratio_range:
        rr = map(float, opts.ratio_range.split(":"))

    run_compare(f, dir_a, dir_b, opts.context, opts.window, opts.binary,
                opts.pvalue, rr)

