import numpy as np
import sys

def moving_window(binfile, window):
    methyl = np.fromfile(binfile, dtype=np.float32)
    mt = np.fromfile(binfile.replace(".methyl.bin",
                                     ".methyltype.bin"),dtype=np.uint8)
    outfile = binfile.replace(".methyl.bin", ".%s.w%i.bin")
    for mname, mtype in (("CG", 1), ("CHG", 2), ("CHH", 3)):
        # make a copy because it gets changed.
        meth = np.copy(methyl) 
        # change the values of the types we're not intrested in to 0
        # the types are 1 - 6 (see README.rst)
        meth[(mt != mtype) & (mt != mtype + 3)] = 0

        assert mt.shape == meth.shape
        of = outfile % (mname, window)
        print "writing: %s" % of
        # only want cg methy
        avg = np.convolve(meth, np.ones((window,)) / float(window),
                      mode='same').astype(np.float32)
        assert avg.shape == mt.shape
        avg.tofile(of)

if __name__ == "__main__":
    usage = """ calculate a moving window average given .methyl.bin files

    usage: %prog -w [window size=50] files

for each input file, 3 output files are created, one for each methylation
type. output files are 32 bit floats each with the same number of entries
as the methyl.bin file (which matches the number of basepairs in the chr).
    """
    import optparse
    p = optparse.OptionParser(usage)
    p.add_option('-w', dest='window', type='int', default=50)
    opts, args = p.parse_args()
    if len(args) == 0:
        print "must specify .methyl.bin files for the moving window"
        sys.exit(p.print_help())

    """
    if len(args) == 1:
        if "*" in args[0]:
            import glob
            args = glob.glob(args[0])
    """
    for binfile in args:
        moving_window(binfile, opts.window)

