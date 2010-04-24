import bsddb
import os

def is_up_to_date_b(a, b):
    if not os.path.exists(b): return False
    return os.stat(a).st_mtime <= os.stat(b).st_mtime

class FileIndex(object):
    ext = ".fidx"

    def is_current(self):
        return is_up_to_date_b(self.filename, self.filename + self.ext)

    @classmethod
    def create(cls, filename, get_next_pos):
        fh = open(filename)
        lines = sum(1 for line in fh)
        bnum = lines if lines > 2**24 else lines * 2
        fh.seek(0)
        db = bsddb.btopen(filename + cls.ext, 'n', cachesize=32768*2,
                          pgsize=512)
        pos = fh.tell()
        while True:
            key = get_next_pos(fh)
            if not key: break
            # always append the | but only used by multiple.
            db[key] = str(pos)
            # fh has been moved forward by get_next.
            pos = fh.tell()
        fh.close()
        db.close()

    def __init__(self, filename, call_class, get_next_pos):
        self.filename = filename
        self.fh = open(self.filename)
        self.call_class = call_class
        if not self.is_current():
            try:
                FileIndex.create(self.filename, get_next_pos)
            except:
                os.unlink(self.filename + self.ext)
                raise

        self.db = bsddb.btopen(filename + self.ext, 'r')

    def __getitem__(self, key):
        # every key has the | appended.
        pos = self.db.get(key, None)
        if pos is None: raise KeyError(key)
        self.fh.seek(long(pos))
        return self.call_class(self.fh)

    def __len__(self):
        return len(self.db)

    def __iter__(self):
        return self.db.iteritems()

    def iterkeys(self):
        return self.db.iterkeys()

if __name__ == "__main__":
    import time
    class FastQEntry(object):
        __slots__ = ('name', 'seq', 'l3', 'qual', 'fpos')
        def __init__(self, fh):
            self.name = fh.readline().rstrip('\r\n')
            self.seq = fh.readline().rstrip('\r\n')
            self.l3 = fh.readline().rstrip('\r\n')
            self.qual = fh.readline().rstrip('\r\n')
 
    f = '/usr/local/src/methylcode/500K.fastq'
    #f = '/usr/local/src/bowtie/bowtie-0.12.1/work/reads/s_1_sequence.txt'

    fi = FileIndex(f, FastQEntry, lambda fh: FastQEntry(fh).name)
    print "items in index:", len(fi)
    N = 10000 
    NKeys = N
    print "getting %i keys..." % N

    it = fi.iterkeys()
    keys = [it.next() for i in xrange(N)]
    import random
    print "shuffling"
    random.shuffle(keys)
    keys = keys[:N]
    print "timing" 
    t = time.time()
    for k in keys:
        entry = fi[k].name
    print time.time() - t
    print len(keys) / (time.time() - t) , "queries per second"
