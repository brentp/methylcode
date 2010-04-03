import sys
import os.path as op
sys.path.insert(0, "/home/brentp/src/methylcode/code/")

from methyl import MethylGroup
prefix = sys.argv[1] # something like: out1234n/thaliana_v9
acontext = sys.argv[2] # CHH or CHG or CG
window = 50

mg = MethylGroup(prefix)


fh = open(mg.dir + mg.prefix + ".test.%ibp.%s.gff" % (window, acontext), "w")
print >>sys.stderr, "writing to %s" % (fh.name, )
print >>fh, "##gff-version 3"
sys.argv[1] = op.abspath(sys.argv[1])
print >>fh, "#%s" % " ".join(sys.argv)

for chr, m in mg.iteritems():
    
    cs, ts, mask = m.as_context(acontext)
    bp_max = len(ts)
    for start in range(0, bp_max + 1, window):
        end = min(start + window, bp_max)
        t_count = ts[start:end].sum()
        c_count = cs[start:end].sum()

        n = mask[start:end].sum()

        if c_count + t_count == 0:
            plot = methyl = 0.0
        else:
            plot = methyl = c_count / float(c_count + t_count)
        strand = "." 
        plot = "%.3g" % plot

        attrs="c=%i;t=%i;n=%i" % (c_count, t_count, n)
        print >>fh, "\t".join(map(str, [chr, sys.argv[0], "dmc", start + 1, end, plot, strand, ".", attrs]))
