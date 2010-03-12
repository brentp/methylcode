from pyfasta import Fasta
import sys
import os
import numpy as np

def check_txt(txt, fa):
    f = Fasta(fa, flatten_inplace=True)
    for line in open(txt):
        seqid, mtype, bp, cs, ts = line.rstrip().split()
        bp = int(bp)
        mtype = int(mtype)
        assert mtype > 0
        if mtype < 4:
            assert f[seqid][bp] == 'C'
        else:
            assert f[seqid][bp] == 'G'
    print txt, "OK"
    return 0

def check_bin(binpath, fa_path):
    # a.5.methyltype.bin => 5
    seqid = binpath.split(".")[-3]
    is_m = ".methyl." in binpath
    is_mtype = ".methyltype." in binpath
    dtype = np.float32 if is_m else np.uint8 if is_mtype else np.uint32

    fa = Fasta(fa_path) 
    bin = np.fromfile(binpath, dtype=dtype)
    assert bin.shape[0] == len(fa[seqid]), (bin.shape[0], len(fa[seqid]))

    assert not np.any(np.isinf(bin))
    assert not np.any(np.isnan(bin))
    assert not np.any(bin < 0)

    if is_m:
        assert 0 == np.min(bin), (binpath, np.min(bin)) 
        assert 1 >= np.max(bin)
        assert 0 < np.average(bin) < 1
    else:
        # TODO: add checks.
        pass
    print binpath, "OK"

if __name__ == "__main__":
    import optparse
    usage = """check output files created by run_bowtie.py
    usage: %prog [options] files_to_check"""
    p = optparse.OptionParser(usage)
    p.add_option("-b", dest="bin", action='store_true',
                 help="check binary files")
    p.add_option("-t", dest="txt", action='store_true',
                 help="check a text file")
    p.add_option("-f", dest="fasta",
                 help="path to the fasta file (required!)")

    opts, args = p.parse_args() 
    if not (opts.bin or opts.txt):
        print "must specify either binary or text file with -b or -t"
        sys.exit(p.print_help())
    if not opts.fasta:
        print "must specify a fasta file"
        sys.exit(p.print_help())

    assert os.path.exists(opts.fasta)


    if opts.txt:
        if not (args[0] and os.path.exists(args[0])):
            print "must specify a txt file to check"
            sys.exit(p.print_help())

        check_txt(args[0], opts.fasta)

    elif opts.bin:
        for binfile in args:
            if not os.path.exists(binfile):
                print "specified binary file: %s does not exist" % binfile
                sys.exit(p.print_help())
            check_bin(binfile, opts.fasta)

