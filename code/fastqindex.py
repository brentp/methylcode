import bsddb
import os

def is_up_to_date_b(a, b):
    if not os.path.exists(b): return False
    return os.stat(a).st_mtime <= os.stat(b).st_mtime

class FastQEntry(object):
    __slots__ = ('name', 'seq', 'l3', 'qual', 'fpos')
    def __init__(self, fh):
        self.name = fh.readline().rstrip('\r\n')
        self.seq = fh.readline().rstrip('\r\n')
        self.l3 = fh.readline()
        self.qual = fh.readline().rstrip('\r\n')

class FastQIndex(object):
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
        while True:
            # fh has been moved forward by get_next.
            pos = fh.tell()
            key = get_next_pos(fh)
            if not key: break
            # always append the | but only used by multiple.
            db[key] = str(pos)
        fh.close()
        db.close()

    def __init__(self, filename, call_class=FastQEntry, get_next_pos=None):
        self.filename = filename
        self.fh = open(self.filename)
        self.call_class = call_class
        if not self.is_current():
            try:
                FastQIndex.create(self.filename, get_next_pos)
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

    def get_position(self, key):
        return long(self.db.get(key, None))

    def __len__(self):
        return len(self.db)

    def __iter__(self):
        return self.db.iteritems()

    def iterkeys(self):
        return self.db.iterkeys()
