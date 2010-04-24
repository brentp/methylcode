"""
given the *same set of reads* aligned to 2 different reference genomes (in the
format of 2 SAM files), find reads that cover a SNP and partition reads into
 2 different sets based on the reference genome for which they match better 
(fewer mismatches).
usage is e.g.:
    python %prog reference_a/methylcoded.sam reference_b/methylcoded.sam

which will write new files:
    reference_a/methylcoded.better.sam
    reference_b/methylcoded.better.sam
which will only contain reads that covered a SNP.

this script requires that tokyocabinet is installed along with the
python module: py-tcdb
"""

import sys
import os.path as op
dir = op.dirname(__file__)
sys.path = [dir, op.join(dir, "lib")] + sys.path
#from utils import SamLine
from fileindex import FileIndex

def fast_sam(fh):
    read = fh.readline()
    fh.readline(); fh.readline(); fh.readline()
    return read.split("\t")[0]

def unique_b(sam_a, sam_b):

    aidx = FileIndex(sam_a, fast_sam, fast_sam)
    out_b = open(sam_b.replace(".sam", ".better.sam"), "w")
    print "writing to %s" % out_b.name
    unique = 0
    for line in open(sam_b):
        b_id = line.split("\t")[0]
        if not b_id in aidx:
            unique += 1
            print >>out_b, line,
    print "unique to %s: %i" % (sam_b, unique)

def main(sam_a, sam_b):
    unique_b(sam_a, sam_b)
    unique_b(sam_b, sam_a)

if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser(__doc__)
    opts, args = p.parse_args()
    if len(args) != 2:
        sys.exit(not p.print_help())

    sam_a, sam_b = args
    main(sam_a, sam_b)
