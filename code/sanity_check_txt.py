from pyfasta import Fasta
def check(txt, fa):
    f = Fasta(fa)
    for line in open(txt):
        seqid, mtype, bp, tot, conv = line.rstrip().split()
        bp = int(bp)
        mtype = int(mtype)
        assert mtype > 0
        if mtype < 4:
            assert f[seqid][bp] == 'C'
        else:
            assert f[seqid][bp] == 'G'
    print "OK"

if __name__ == "__main__":
    import sys
    txt = sys.argv[1]
    fa  = sys.argv[2]

    check(txt, fa)
