import numpy as np
cimport numpy as np
cimport stdlib
import sys


def _update_conversions(char *ref_seq, char *aln_seq, int base_position, 
                        char *pair,
                        np.ndarray[np.uint32_t] c_count,
                        np.ndarray[np.uint32_t] t_count,
                       int allowed_mismatches, int n):
    """
    this updates the conversions (counts) from c2t in between
    the ref_sequence and the aln_seq
    first it checks if the sequence was mistakenly aligned by the c2t 
    conversion. if so, it does not update the counts.
    """
    cdef int i #, n = stdlib.strlen(ref_seq)
    cdef char a, b
    # CT, GA
    cdef char c1 = pair[0], c2 = pair[1]
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
    d = ["."] * n
    """

    for i in range(n):
        a = ref_seq[i]
        if a != c1: continue
        b = aln_seq[i]
        # CC
        if b == c1:
            c_count[base_position + i] = 1
            # d[i] = chr(c1)
        # CT
        if b == c2:
            t_count[base_position + i] = 1
            # d[i] = chr(c2)
    # sys.stderr.write("mat:" + "".join(d) + "\n")
    return 0
