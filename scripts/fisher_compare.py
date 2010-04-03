from fisher import pvalue_npy, pvalue
import sys
import os.path as op
sys.path.insert(0, op.join(op.dirname(__file__), "../code"))
from methyl import MethylGroup
import numpy as np

def bin_setup(chr, adir, bdir, context='CG'):
    am = MethylGroup(adir)[chr]
    bm = MethylGroup(bdir)[chr]
    return am.as_context(context), bm.as_context(context)

def difference_compare(w_ts, w_cs, m_ts, m_cs):
    w_meth = w_cs.astype('f')/ (w_cs + w_ts)
    m_meth = m_cs.astype('f')/ (m_cs + m_ts)
    m_diff = w_meth - m_meth
    m_diff[np.isnan(m_diff)] = 0.0

    window_size = 100.0
    print "convolving..."
    avg_win = np.convolve(m_diff, np.ones(window_size, 
                                 dtype='f')/window_size, mode="same")
    print "sorting..."
    b = np.sort(avg_win)
    cutoff = b[int(0.995 * b.shape[0])]
    print cutoff
    min_idx = np.where(b >= cutoff)[0][0]
    different_idxs, = np.where(avg_win > b[min_idx])

    first = different_idxs[0]
    previous = first
    for last in different_idxs:
        if last - previous > window_size:
            yield first, previous
            first = last
        previous = last

    yield first, previous

def fisher_compare(w_ts, w_cs, m_ts, m_cs):
    left_tail, right_tail, two_tail = pvalue_npy(w_cs.astype(np.uint), 
                                                 w_ts.astype(np.uint),  
                                                 m_cs.astype(np.uint), 
                                                 m_ts.astype(np.uint))

    tail = right_tail
    p05 = tail < 0.05

    window_size = 100.0
    avg_win = np.convolve(p05, np.ones(window_size, 
                                 dtype='f')/window_size, mode="same")

    b = np.sort(avg_win)
    # so next lines find the index in b of the *first* value that
    # is at/above the 99.5th pctile.
    cutoff = b[int(0.995 * b.shape[0])]
    min_idx = np.where(b >= cutoff)[0][0]
    different_idxs, = np.where(avg_win > b[min_idx])
    """ 
    # NOTE TEMPORARY CHANGE
    """
    different_idxs, = np.where(p05)
    for idx in different_idxs + 1:
        yield idx, tail[idx - 1]

    """
    first = different_idxs[0]
    previous = first
    for last in different_idxs:
        if last - previous > window_size:
            yield first, previous
            first = last
        previous = last

    yield first, previous
    """

"""
1   TAIR9   gene    3631    5899    .   +   .   ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
1   TAIR9   mRNA    3631    5899    .   +   .   ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
"""


def run_gff(flat, mpatt, wpatt):
    fh = open('fisher.different.bp.gff', 'w')
    for i in range(1, 6):
        chr = str(i)
        w_ts, w_cs, m_ts, m_cs = bin_setup(chr, mpatt, wpatt)
        print >>fh, "##gff-version 3"
        for start1, pv in fisher_compare(w_ts, w_cs, m_ts, m_cs):
            strand = flat.fasta[chr][start1 - 1]
            assert strand in "GC", strand
            strand = "+" if strand == "C" else "-"
            attrs="pvalue=%.3g" % pv
            accns = flat.get_features_in_region(chr, start1, start1 + 1)
            accns = [a["accn"] for a in accns]
            if accns:
                attrs +=";accns=" + ",".join(accns)
            print >>fh, "\t".join([chr, "fischerlab", "dmc", str(start1), str(start1), ".", strand, ".", attrs])

def run_accn_gff(flat, mpatt, wpatt):
    fh = open('fisher.different.accn.gff', 'w')
    print >>fh, "##gff-version 3"
    for i in range(1, 6):
        chr = str(i)
        f = flat[flat['seqid'] == chr]
        w_ts, w_cs, m_ts, m_cs = bin_setup(chr, mpatt, wpatt)
        for row in f:
            slicer = slice(row['start'] - 1, row['end'])
            w_t_count = w_ts[slicer].sum()
            w_c_count = w_cs[slicer].sum()
            m_t_count = m_ts[slicer].sum()
            m_c_count = m_cs[slicer].sum()

            p = pvalue(w_t_count, w_c_count, m_t_count, m_c_count)
            pv = p.two_tail
            if pv > 0.01: continue
            gc = f.fasta[chr][row['start'] - 1: row['end']].upper()
            gc = gc.count("G") + gc.count("C")
            attrs="ID=%s;p=%.3g;wtc=%i;wtt=%i;mtc=%i;mtt=%i;gc=%i" % \
                        (row['accn'], pv, w_t_count, w_c_count, m_t_count, m_c_count, gc)
            print >>fh, "\t".join(map(str, [chr, "fischerlab", "gene", row['start'], row['end'], ".", row['strand'], ".", attrs]))

def run_50bp_gff(flat, adir, bdir, context, window=50):
    fh = open('fisher.different.%s.%ibp.gff' % (context, window), 'w')
    print >>sys.stderr, "writing to:", fh.name
    print >>fh, "##gff-version 3"
    for i in range(1, 6):
        chr = str(i)
        bp_max = len(flat.fasta[chr])
        (a_cs, a_ts, a_mask), (b_cs, b_ts, b_mask) = bin_setup(chr, adir, bdir, context)
        for start in range(0, bp_max + window, window):
            end = min(start + window, bp_max)
            a_t_count = a_ts[start:end].sum()
            a_c_count = a_cs[start:end].sum()
            b_t_count = b_ts[start:end].sum()
            b_c_count = b_cs[start:end].sum()

            p = pvalue(a_t_count, a_c_count, b_t_count, b_c_count)
            pv = p.two_tail
            if pv > 0.01: continue
            gc = f.fasta[chr][start:end].upper()
            gc = gc.count("G") + gc.count("C")

            a_methyl = a_c_count / float(a_c_count + a_t_count)
            b_methyl = b_c_count / float(b_c_count + b_t_count)
            strand = "+" if a_methyl > b_methyl else "-"
            strand = "."
            plot = a_methyl - b_methyl
            attrs="p=%.3g;ac=%i;at=%i;bc=%i;bt=%i;gc=%i" % \
                        (pv, a_c_count, a_t_count, b_c_count, b_t_count, gc)
            accns = flat.get_features_in_region(chr, start + 1, end)
            accns = [a["accn"] for a in accns]
            if accns:
                attrs +=";accns=" + ",".join(accns)
            print >>fh, "\t".join(map(str, [chr, "methylation", "dmc", start + 1, end, plot, strand, ".", attrs]))
    

def run_all(flat, mpatt, wpatt):

    all_f = []
    all_d = []
    fh_f = open('fisher.different.txt', 'w')
    fh_d = open('difference.different.txt', 'w')
    for i in range(1, 6):
        chr = str(i)
        w_ts, w_cs, m_ts, m_cs = bin_setup(chr, mpatt, wpatt)

        for start0, end0 in fisher_compare(w_ts, w_cs, m_ts, m_cs):
            accns = flat.get_features_in_region(chr, start0 + 1, end0 + 1)
            accns = [a["accn"] for a in accns]
            all_f.extend(accns) 
            print >>fh_f, chr, start0, end0, ",".join(accns)

        for start0, end0 in difference_compare(w_ts, w_cs, m_ts, m_cs):
            accns = flat.get_features_in_region(chr, start0 + 1, end0 + 1)
            accns = [a["accn"] for a in accns]
            all_d.extend(accns) 
            print >>fh_d, chr, start0, end0, ",".join(accns)

    print >>open('fisher.different.accns', 'w'), "\n".join(sorted(list(set(all_f))))
    print >>open('difference.different.accns', 'w'), "\n".join(sorted(list(set(all_d))))

if __name__ == "__main__":
    from flatfeature import Flat
    #f = Flat("/opt/src/athly/data/thaliana_v9.flat", "/opt/src/athly/data/thaliana_v9.fasta")
    f = Flat("/labdata/thaliana_v9/thaliana_v9.flat", "/labdata/thaliana_v9/thaliana_v9.fasta")
    import sys
    #run_all(f, "out5678n/thaliana_v9.%s.converted.bin",
    #           "out1234n/thaliana_v9.%s.converted.bin")
    #run_gff(f, "out5678n/thaliana_v9.%s.converted.bin",
    #           "out1234n/thaliana_v9.%s.converted.bin")
    #run_accn_gff(f, "out5678n/thaliana_v9.%s.converted.bin",
    #                "out1234n/thaliana_v9.%s.converted.bin")
    if len(sys.argv) > 4:
        window = int(sys.argv[4])
    else:
        window = 50
    run_50bp_gff(f, sys.argv[1], sys.argv[2], sys.argv[3].upper(), window)

