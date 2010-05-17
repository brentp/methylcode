"""
check a fastq file for validity checks are:
   1) 4 line rrcords
   2) header starts with @
   3) qual line-length == seq line-length
"""

import sys

fh = open(sys.argv[1])
rl = fh.readline
while True:
    header = rl().strip()
    if not header: continue
    assert header[0] == "@", (header)
    seq = rl()
    assert rl()[0] == "+"
    qual = rl()
    assert len(seq) == len(qual), (header, seq, qual)
