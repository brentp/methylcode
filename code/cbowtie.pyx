import numpy as np
cimport numpy as np
cimport stdlib
import sys


def _update_conversions(char *ref_seq, char *aln_seq, int base_position, 
                        char *pair,
                        np.ndarray[np.uint8_t] total_count, 
                        np.ndarray[np.uint8_t] converted_count,
                       int allowed_mismatches):
    """
    this updates the conversions (counts) from c2t in between
    the ref_sequence and the aln_seq
    first it checks if the sequence was mistakenly aligned by the c2t 
    conversion. if so, it does not update the counts.
    """
    cdef int i, n = stdlib.strlen(ref_seq)
    cdef char a, b
    # CT, GA
    cdef char c1 = pair[0], c2 = pair[1]

    """
    for i in range(n):
        # TODO: check this for GA
        if ref_seq[i] != c2: continue
        if aln_seq[i] == c1: 
            if allowed_mismatches == 0: return 1
            allowed_mismatches -= 1
    """

    for i in range(n):
        a = ref_seq[i]
        if a != c1: continue
        b = aln_seq[i]
        total_count[base_position + i] += 1
        if b == c2:
            converted_count[base_position + i] += 1
    return 0
