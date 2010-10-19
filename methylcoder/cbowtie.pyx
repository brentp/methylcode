#cython: boundscheck=False
#cython: wraparound=False

from libc.stdlib cimport malloc, free
from libc.string cimport strlen

import numpy as np
cimport numpy as np
import sys


def _update_conversions(char *ref_seq, char *aln_seq, int base_position,
                        char *pair,
                        np.ndarray[np.uint32_t] c_count,
                        np.ndarray[np.uint32_t] t_count,
                       int allowed_mismatches, int n,
                       bint is_colorspace):
    """
    this updates the conversions (counts) from c2t in between
    the ref_sequence and the aln_seq
    first it checks if the sequence was mistakenly aligned by the c2t 
    conversion. if so, it does not update the counts.
    """
    cdef int i
    cdef char a, b
    if is_colorspace:
        t = cs2seq(aln_seq)
        aln_seq = t
    # CT, GA
    cdef char c1 = pair[0], c2 = pair[1]
    if allowed_mismatches > -1:
        for i in range(n):
            # if genome is a T and read is a C
            # bowtie couldnt find this cause we
            # converted C to T.
            if ref_seq[i] != c2: continue
            if aln_seq[i] == c1:
                if allowed_mismatches == 0: return 1
                allowed_mismatches -= 1
        if allowed_mismatches < 0: return 1

    # for proofing.
    """
    sys.stderr.write("ref:" + ref_seq + "\n")
    sys.stderr.write("aln:" + aln_seq + "\n")
    tts = 0
    ccs = 0
    d = ["."] * n
    """

    for i in range(n):
        a = ref_seq[i]
        if a != c1: continue
        b = aln_seq[i]
        # CC
        if b == c1:
            c_count[base_position + i] += 1
            # debug
            #d[i] = chr(c1); ccs += 1
        # CT
        elif b == c2:
            t_count[base_position + i] += 1
            # debug
            #d[i] = chr(c2); tts += 1
    """
    sys.stderr.write("mat:" + "".join(d) + "\n")
    sys.stderr.write("remained     c: %i\n" % ccs)
    sys.stderr.write("converted to t: %i\n\n" % tts)
    """
    return 0

cdef dict colorspace = {
    65: [ # A
        65, # 0
        67, # 1
        71, # 2
        84, # 3
        78, # 4
        78],# .
    67: # C
        [ 67, 65, 84, 71, 78, 78],
    71: # G
        [71, 84, 65, 67, 78, 78],
    84:  # T
        [84, 71, 67, 65, 78, 78],
    78:  # N
        [ 78, 78, 78, 78, 78],
}


def cs2seq(cs):
    """
    convert a color-space encoded read to DNA sequence
    >>> cs2seq('A32102302224')
    'ATCAAGCCTCTN'

    >>> cs2seq('42220320123A')
    'NTCTCCGAACTA'

    """
    if not cs[0] in "ACGTN":
        cs = cs[::-1]
        return _cs2seq(cs)[::-1]
    return _cs2seq(cs)


cdef _cs2seq(char *cs):
    cdef int seqlen = strlen(cs)
    #print >>sys.stderr, "\n\n%s" % cs
    cdef int i = 0
    cdef char *seq = <char *>malloc(sizeof(char) * (1 + seqlen))
    seq[0] = cs[0]
    #cdef object k1, k2
    cdef char k1, k2
    cdef int start = 1, stop = seqlen
    cdef list sublist

    for i in range(start, stop):

        k1 = seq[i - 1] # decoded DNA character
        k2 = cs[i] - 48 # char number to index
        sublist = colorspace[k1]
        seq[i] = sublist[k2]
    seq[seqlen] = c'\0'
    try:
        return seq
    finally:
        free(seq)
