import sys

class FastQ(object):
    def __init__(self, fname, unique=True, max_n=8):
        self.fname = fname
        self.uniq = unique
        self.max_n = max_n

    def __iter__(self):
        self.fh = open(self.fname)
        self.seen = {}
        return self

    def next(self):

        while True:
            name = self.fh.readline().rstrip()
            if not name: raise StopIteration
            assert name[0] == "@"

            fpos = self.fh.tell()
            seq = self.fh.readline().rstrip()
            l3 = self.fh.readline()[0] # redundant, just use single char.
            qual = self.fh.readline().rstrip()
            if self.uniq:
                if name in self.seen: continue
                self.seen[name] = None
            if self.max_n:
                if seq.count("N") > self.max_n: continue

            return FastQEntry(name, seq, l3, qual, fpos)

class FastQEntry(object):
    __slots__ = ('name', 'seq', 'l3', 'qual', 'fpos')
    def __init__(self, name, seq, l3, qual, fpos):
        self.fpos = fpos
        self.name = name
        self.seq = seq
        self.l3 = l3
        self.qual = qual

    def __str__(self):
        return "\n".join((self.name, self.seq, self.l3, self.qual))
    
    def __repr__(self):
        return "F(%s: %i)" % (self.name, self.fpos)


if __name__ == "__main__":
    import optparse
    p = optparse.OptionParser("%prog [options] fastq")
    p.add_option("-n", dest="max_n", type='int', default=8, 
                 help="filter lines with more N's than this")
    p.add_option("-u", dest="unique", action="store_true", default=False, 
                 help="only print unique sequences")
    opts, args = p.parse_args()
    if not len(args):
        print "send in fastq file!"
        sys.exit(p.print_help())

    for q in FastQ(args[0], opts.unique, opts.max_n):
        print q
