#cython: boundscheck=False
#cython: wraparound=False

from libc.stdlib cimport malloc, free, atoi
from libc.string cimport strlen, strchr, strncpy

import numpy as np
cimport numpy as np
import sys

cdef inline void set_soft_masking_start_stop(char *cigar, int n, int *s1i, int *s2i):
    cdef char *M, *S1, *S2
    cdef char[3] tmp_s
    if strchr(cigar, c'S') == NULL: return

    M = strchr(cigar, c'M')
    S2 = strchr(M, c'S')
    S1 = strchr(cigar, c'S')
    if S1 == S2:
        S1 = NULL
        s1i[0] = 0
    else:
        tmp_s[2] = tmp_s[1] = c'\0'
        strncpy(tmp_s, cigar, S1 - cigar)
        s1i[0] = atoi(tmp_s)

    tmp_s[1] = c'\0'
    strncpy(tmp_s, M + 1, 2)
    tmp_s[2] = '\0'
    # it's stuff skipped at the end. so we keep up to there.
    s2i[0] = (n - atoi(tmp_s))
    #assert *s1i < *s2i, (cigar, M + 1, S1, S2, tmp_s, atoi(tmp_s))



def _update_conversions(char *ref_seq, char *aln_seq, int base_position,
                        char *pair,
                        np.ndarray[np.uint32_t] c_count,
                        np.ndarray[np.uint32_t] t_count,
                       int allowed_mismatches, int n,
                       char *cigar):
    """
    this updates the conversions (counts) from c2t in between
    the ref_sequence and the aln_seq
    first it checks if the sequence was mistakenly aligned by the c2t 
    conversion. if so, it does not update the counts.
    use cigar to account for the 'S'oft masking.
    """
    cdef int i
    cdef char a, b
    cdef char c1 = pair[0], c2 = pair[1]

    ######################################################
    # handle soft masking in cigar
    ######################################################
    cdef int s1i = 0, s2i = n
    set_soft_masking_start_stop(cigar, n, &s1i, &s2i)
    ######################################################

    if allowed_mismatches > -1:
        for i in range(s1i, s2i):
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
    sys.stderr.write("cigar:" + cigar + "\n")
    sys.stderr.write("ref:" + (" " * s1i) + ref_seq + "\n")
    sys.stderr.write("aln:" + aln_seq + "\n")
    tts = 0
    ccs = 0
    d = ["."] * n
    for i in range(s1i):
        d[i] = " "
    """

    for i in range(0, s2i):
        if i < s1i: continue
        a = ref_seq[i - s1i]
        if a != c1: continue
        b = aln_seq[i]
        # CC
        # TODO: check these base positions!!!
        if b == c1:
            c_count[base_position + i - s1i] += 1
            # debug d[i] = chr(c1); ccs += 1
        # CT
        elif b == c2:
            t_count[base_position + i - s1i] += 1
            # debug d[i] = chr(c2); tts += 1
    """
    for i in range(s2i, n):
        d[i] = " "
    sys.stderr.write("mat:" + "".join(d) + "\n")
    sys.stderr.write("remained     c: %i\n" % ccs)
    sys.stderr.write("converted to t: %i\n\n" % tts)
    raw_input("... press a key ...")
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

cdef dict seqspace = {
    'AT': 3,
    'AG': 2,
    'AC': 1,
    'AA': 0,
    'CT': 2,
    'CG': 3,
    'CC': 0,
    'CA': 1,
    'GT': 1,
    'GG': 0,
    'GC': 3,
    'GA': 2,
    'TT': 0,
    'TG': 1,
    'TC': 2,
    'TA': 3,
    }

def seq2cs(char *seq):
    """
    convert sequence into colorspace
    """
    cdef int seqlen = strlen(seq)
    cdef char *cs = <char *>malloc(sizeof(char) * (1 + seqlen))
    cdef char[3] dub
    dub[2] = '\0'
    cs[0] = seq[0]
    cdef int i
    for i in range(seqlen - 1):
        if seq[i] == 78: cs[i + 1] = '4'
        elif seq[i + 1] == 78: cs[i + 1] = '4'
        else:
            dub[0] = seq[i]
            dub[1] = seq[i + 1]
            cs[i + 1] = seqspace[dub] + 48
    cs[seqlen] = c'\0'
    try:
        return cs
    finally:
        free(cs)

