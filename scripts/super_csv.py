"""
given the binary files created with run_bowtie.py and a .gff
create a super csv by gene with 136 columns. (goes to stdout when this script
is run. adjust the paths at the bottom before running. requires 
flatfeature/fatfeature.py in bpbio repository

the 136 is broken down in to 1 + (3 * 5 * 9) groups
where (i use "context" to mean cg, chg, or chh):
    1 is the accn name
    3 is the 3 types of methylation context
    5 is for the methylation stats that are reported [see below]
    9 is for the different locs (gene, cds, intron, upstream 10/100/100,
                                and downstream 10/100/100) 
    (up and downstream exclude any intervening CDSs).

the 5 stats reported are:
    + avg : average methylation for sites in the given context (cg or chg or chh)
    + avg_gt0 : average methylaiton for sites in this context which have
                any methylation (exclude 0's--this seems to be a common stat used in
                the methylation papers)
    + n_methylated: the number of sites in this context that can be
                    methylated and are
    + n: the number of sites in this context that can be methylated.
    + %_methylated: the ratio of the 2 previous columns == the proportion
                    of sites that can be methylated in this context that are. so it's
                    actually proportion, not %.

so the 136 headers have the format:

[type]_[context]_[stat]
eg:
[gene]_[cg]_[%_methylated] => gene_cg_%_methylated
"""

from fatfeature import Fat
import numpy as np
import sys
import os


def pairs_to_slice(pairs):
    """
    given a list of tuples (like a list of CDS start, stops), return
    the numpy array that will work as a slice for those tuples
    """
    return np.concatenate([np.arange(s0-1, s1) for s0, s1 in pairs])


def calc_stats(mfiles_pat, gff):
    fat = Fat(gff)

    methyl = {}
    contexts = {}
    for i in range(1, 6):
        i = str(i)
        mf = mfiles_pat % i
        assert os.path.exists(mf)
        methyl[i] = np.fromfile(mf, dtype=np.float32)
        mt = np.fromfile(mf.replace(".methyl.", ".methyltype."), dtype=np.uint8)
        # these can be used to mask to a given context.
        contexts[i] = {'cg': (mt == 1) | (mt == 4),
                       'chg': (mt == 2) | (mt == 5),
                       'chh': (mt == 3) | (mt == 6)}
    header = ["accn"]
    for ctx in ('cg', 'chg', 'chh'):
        # TODO: make this suck less.
        header.append(",".join([
            "gene_CTX_avg,gene_CTX_avg_gt0,gene_CTX_n_methylated,gene_CTX_n,gene_CTX_%_methylated",
            "cds_CTX_avg,cds_CTX_avg_gt0,cds_CTX_n_methylated,cds_CTX_n,cds_CTX_%_methylated",
            "intron_CTX_avg,intron_CTX_avg_gt0,intron_CTX_n_methylated,intron_CTX_n,intron_CTX_%_methylated",
            "up10_CTX_avg,up10_CTX_avg_gt0,up10_CTX_n_methylated,up10_CTX_n,up10_CTX_%_methylated",
            "up100_CTX_avg,up100_CTX_avg_gt0,up100_CTX_n_methylated,up100_CTX_n,up100_CTX_%_methylated",
            "up1000_CTX_avg,up1000_CTX_avg_gt0,up1000_CTX_n_methylated,up1000_CTX_n,up1000_CTX_%_methylated",
            "down10_CTX_avg,down10_CTX_avg_gt0,down10_CTX_n_methylated,down10_CTX_n,down10_CTX_%_methylated",
            "down100_CTX_avg,down100_CTX_avg_gt0,down100_CTX_n_methylated,down100_CTX_n,down100_CTX_%_methylated",
            "down1000_CTX_avg,down1000_CTX_avg_gt0,down1000_CTX_n_methylated,down1000_CTX_n,down1000_CTX_%_methylated"
        ]).replace('CTX', ctx)) # bleckh. shrug.
    print ",".join(header)

    for accn, f in sorted(fat.iteritems()):
        
        data = [accn]
        if not f.seqid in contexts: continue # C, G

        for mtype, context in sorted(contexts[f.seqid].iteritems()):
            for locs in ([[f.start, f.end]], 
                         getattr(f, 'CDS', None), 
                         fat.introns(f),
                         fat.upstream(f, 10, noncoding=True), 
                         fat.upstream(f, 100, noncoding=True), 
                         fat.upstream(f, 1000, noncoding=True),
                         fat.downstream(f, 10, noncoding=True), 
                         fat.downstream(f, 100, noncoding=True), 
                         fat.downstream(f, 1000, noncoding=True)
            ):
                if locs is None:
                    # occurs when there's no CDS.
                    data.extend(["na","na","na","na","na"])
                    continue

                slicer = pairs_to_slice(locs)
                try:
                    ctx = context[slicer] # this context.
                except IndexError: 
                    # difference between fasta and features due to version
                    slicer = slicer[slicer < context.shape[0]]
                    ctx = context[slicer]

                # methylation for this CDS masked to current context
                me = methyl[f.seqid][slicer] * ctx 
                # number of sites in this context.
                ctx_sites = ctx.sum() 
                # number of site in this context that are methylated.
                m_ctx_sites = (me > 0).sum() 
                # proprtion of sites that can be methylated that are.
                p_methylated = float(m_ctx_sites) / ctx_sites

                # average methylation for sites in this context.
                avg_methyl = me.sum() / float(ctx_sites)

                # average methylation for sites that are methylated. (exclude zeros).
                avg_methyl_gt0 = me.sum() / float(m_ctx_sites)
                data.extend(["%.5f" % d for d in [avg_methyl, avg_methyl_gt0]])
                data.extend(["%i" % d for d in [m_ctx_sites, ctx_sites]])
                data.append("%.5f" % p_methylated)
        
        print ",".join(data)



if __name__ == "__main__":

    # adjust this accordingly.
    methpat = "/usr/local/src/methylcode/out1234n/thaliana_v9.%s.methyl.bin"

    gff = "/labdata/thaliana_v9/thaliana_v9_genes.gff"
    calc_stats(methpat, gff)
