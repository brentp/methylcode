"""
given a fasta file. for each chromosome, write a file containing:
    0: no methylation possible (not a C or G)
    1: + CG methylation 
    2: + CHG methylation
    3: + CHH methylation
    4: - CG methylation 
    5: - CHG methylation
    6: - CHH methylation
"""
import numpy as np
import string

_complement = string.maketrans('ATCG', 'TAGC')
revcomp = lambda s: s.translate(_complement)[::-1]

# + 1   2    3
[_, CG, CHG, CHH] = range(4)
# - 4   5    6

def calc_methylation(sequence):
    """
    >>> calc_methylation("CCC")
    array([3, 3, 3], dtype=uint8)

    >>> calc_methylation("CCCCC")
    array([3, 3, 3, 3, 3], dtype=uint8)

    >>> calc_methylation("GGGGG")
    array([6, 6, 6, 6, 6], dtype=uint8)

    >>> calc_methylation("GGCGG")
    array([6, 6, 1, 4, 5], dtype=uint8)

    >>> calc_methylation("GGAGG")
    array([6, 6, 0, 6, 6], dtype=uint8)
    """
    sequence = "HH" + sequence.upper() + "HH"

    seq = np.array(sequence, dtype='c')
    methyl = _calc_methylation(seq)

    seq[:] = np.array(revcomp(sequence), dtype='c')
    methyl = _calc_methylation(seq, methyl[::-1])
    # chop off the 'HH' and reverse since it was reversed
    # for the 2nd run. the index is confusing, but correct.
    return methyl[-3:1:-1]


def _calc_methylation(seq, methyl=None):
    """
    NOTE: this does not handle the case where a 'C'
    appears within 2bp of the end (it will try to index)
    above when checking if subsequent basepairs are a  G
    that correction is handled in calc_methylation

    >>> _calc_methylation(np.array("CAG", dtype='c'))
    array([2, 0, 0], dtype=uint8)

    # and the reverse complement:
    >>> _calc_methylation(np.array("GTCCGG", dtype='c'))
    array([0, 0, 2, 1, 0, 0], dtype=uint8)

    """
    # to differentiate between plus and minus
    # methylation. plus is CG,CHG, CHH => 1, 2, 3
    # and minus is 4, 5, 6
    # adder is 0 for + strand and 3 for -
    adder = 3
    if methyl is None:
        methyl = np.zeros(seq.shape, dtype=np.uint8)
        adder = 0

    # CG
    c_idxs, = np.where(seq == 'C')
    # where a 'G' follows a 'C'
    first_g = seq[c_idxs + 1] == 'G'
    cg_idxs, = np.where(first_g)
    cg = c_idxs[cg_idxs]
    methyl[cg] = CG + adder
    del cg

    # CHG
    second_g = seq[c_idxs + 2] == 'G'
    chg_idxs, = np.where((~first_g) & (second_g))
    chg = c_idxs[chg_idxs]
    methyl[chg] = CHG + adder
    del chg

    # CHH
    chh_idxs, = np.where(~(first_g | second_g))
    chh = c_idxs[chh_idxs]
    methyl[chh] = CHH + adder

    return methyl

def write_file(seq, file_path):
    methyl = calc_methylation(str(seq))
    methyl.tofile(file_path)


if __name__ == "__main__":
    import optparse
    import sys
    import doctest
    import os
    doctest.testmod()
    import sys
    sys.exit()

    p = optparse.OptionParser(__doc__)
    opts, args = p.parse_args()
    if not (len(args) and args[0]): sys.exit(p.print_help())

    from pyfasta import Fasta
    path = args[0]
    f = Fasta(path)
    for k in f.iterkeys():
        print "calculating methylation data for seqid: %s" % k
        binpath = "%s.%s.methyl.type.bin" % (path, k)
        write_file(f[k], binpath)
        print "saved to %s" % binpath
