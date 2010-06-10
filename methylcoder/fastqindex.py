import bsddb
import os

def is_up_to_date_b(a, b):
    if not os.path.exists(b): return False
    return os.stat(a).st_mtime <= os.stat(b).st_mtime

class FastQEntry(object):
    lines = 4
    bowtie_flag = 'q'
    __slots__ = ('name', 'seq', 'l3', 'qual', 'fpos')
    def __init__(self, fh):
        self.name = fh.readline().rstrip('\r\n')
        self.seq = fh.readline().rstrip('\r\n')
        self.l3 = fh.readline()
        self.qual = fh.readline().rstrip('\r\n')

    @classmethod
    def get_next_position(kls, fh):
        key = fh.readline().rstrip("\n")
        for i in range(kls.lines - 1): fh.readline()
        return key

def guess_index_class(filename):
    assert os.path.exists(filename)
    fh = open(filename)
    header = fh.readline()
    assert header[0] in '@>', (header, "should be fastq or fasta")
    if header[0] == "@": return FastQIndex
    return FastaIndex


class FastaEntry(FastQEntry):
    """
    NOTE: this only handles 1 line sequence in the fasta file.
    """
    lines = 2
    bowtie_flag = 'f'
    __slots__ = ('name', 'seq')
    def __init__(self, fh):
        self.name = fh.readline().rstrip('\r\n')
        self.seq = fh.readline().rstrip('\r\n')

class RawEntry(FastQEntry):
    lines = 1
    bowtie_flag = 'r'
    __slots__ = ('name', 'seq')
    def __init__(self, fh):
        self.seq = fh.readline().rstrip('\r\n')


class FastQIndex(object):
    ext = ".fqdx"
    entry_class = FastQEntry

    def is_current(self):
        return is_up_to_date_b(self.filename, self.filename + self.ext)

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
            db[key] = str(pos)
        fh.close()
        db.close()

    def __init__(self, filename):
        self.filename = filename
        self.fh = open(self.filename)
        if not self.is_current():
            try:
                FastQIndex.create(self.filename, self.entry_class.get_next_position)
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

class FastaIndex(FastQIndex):
    ext = '.fadx'
    entry_class = FastaEntry
    def __init__(self, filename):
        FastQIndex. __init__(self, filename)
