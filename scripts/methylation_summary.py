import sys
import os.path as op
sys.path.insert(0, op.join(op.dirname(__file__), "../code"))
from methyl import MethylGroup
mg = MethylGroup(sys.argv[1])
contexts = ('CG', 'CHG', 'CHH')

 
print "#", mg.prefix, mg.dir, mg.pattern
 
total_cs = dict((ctx, 0) for ctx in contexts)
total_ts = dict((ctx, 0) for ctx in contexts)

M = dict((ctx, 0) for ctx in contexts)
P = dict((ctx, 0) for ctx in contexts)

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
        n_methylated = (cg_cs > 0).sum()
        possible_methylated = cg_mask.sum()
        P[context] += possible_methylated
        M[context] += n_methylated
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


ymax=M['CHH'] + P['CHH']
CGPCT=int(100 * M['CG'] / P['CG'] + 0.5)
CHGPCT=int(100 * M['CHG'] / P['CHG'] + 0.5)
CHHPCT=int(100 * M['CHH'] / P['CHH'] + 0.5)
chart = """\
cht=bvs
chs=350x350
chbh=a
chco=ff0000|00ff00|0000ff,ff000088|00ff0088|0000ff88
chd=t:%(CGM)i,%(CHGM)i,%(CHHM)i|%(CGP)i,%(CHGP)i,%(CHHP)i
chxt=x,y
chxl=0:|CG|CHG|CHH
chxr=1,0,%(ymax)i
chds=0,%(ymax)i
chtt=methlated+/+total+sites
chm=tMethylated,000000,0,2,12|tTotal+Sites,000000,1,2,12
|t%(CGPCT)i%%,000000,1,0,12,1,c::1
|t%(CHGPCT)i%%,000000,1,1,12,1,c::1
|t%(CHHPCT)i%%,000000,1,2,12,1,c::1
""" % dict(CGM=M['CG'], CHGM=M['CHG'], CHHM=M['CHH'], \
           CGP=P['CG'], CHGP=P['CHG'], CHHP=P['CHH'],
           ymax=ymax, CGPCT=CGPCT, CHGPCT=CHGPCT, CHHPCT=CHHPCT)

chart = 'http://chart.apis.google.com/chart?' + chart.replace("\n|", "").replace("\n", "&")
print chart
