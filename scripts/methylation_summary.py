import sys
import os.path as op
sys.path.insert(0, op.join(op.dirname(__file__), "../code"))
from methyl import MethylGroup
mg = MethylGroup(sys.argv[1])
contexts = ('CG', 'CHG', 'CHH')
 
print "#", mg.prefix, mg.dir, mg.pattern
 
total_cs = dict((ctx, 0) for ctx in contexts)
total_ts = dict((ctx, 0) for ctx in contexts)
print "seqid,context,p_methylated,total_possible_sites,possible_sites_covered_by_reads,cs,ts,cs/(cs + ts)"
for seqid, meth in mg.iteritems():
    for context in contexts:
        if seqid in 'CcMm': continue
        cg_cs, cg_ts, cg_mask = meth.as_context(context)
        total_sites = cg_mask.sum()
        mask = (cg_cs + cg_ts) > 0
        cg_cs = cg_cs[mask]
        cg_ts = cg_ts[mask]
        cg_mask = cg_mask[mask]
        methylation = cg_cs.astype('f') / (cg_ts + cg_cs)
        n_methylated = (methylation > 0.0).sum()
        possible_methylated = cg_mask.sum()
        proportion_methylated = float(n_methylated) / possible_methylated
        rat = float(cg_cs.sum())
        rat /= (rat + cg_ts.sum())
        total_ts[context] += cg_ts.sum()
        total_cs[context] += cg_cs.sum()
        rat = "%.5f" % rat
        proportion_methylated = "%.5f" % proportion_methylated

        print ",".join(map(str, (seqid, context, proportion_methylated, 
                     total_sites, cg_mask.sum(), cg_cs.sum(), cg_ts.sum(), rat)))

print '\n# genome wide methylation (c /(c + t))'
for context in contexts:
    cs = total_cs[context]
    ts = total_ts[context]
    print '# %s: %.4f'  % (context, cs / float(cs + ts))
