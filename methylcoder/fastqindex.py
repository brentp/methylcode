import os
try:
    import bsddb
    def get_db(dbname, mode):
        return bsddb.btopen(dbname, mode, cachesize=32768*2,
                          pgsize=512)

except ImportError:
    from tcdb.bdb import BDBSimple as BDB
    import tcdb.bdb as tc

    def get_db(dbname, mode):
        db = BDB()
        if mode in 'cn':
            db.open(dbname,
                    omode=tc.OWRITER | tc.OTRUNC | tc.OCREAT,
                    bnum=int(1e7), lcnum=2**19,
                    apow=6, opts=tc.TLARGE, xmsiz=2**26)
        else:
            db.open(dbname, omode=tc.OREADER)
        return db

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

    def __str__(self):
        return "\n".join((self.name, self.seq, self.l3, self.qual))

def guess_index_class(filename):
    assert os.path.exists(filename)
    fh = open(filename)
    advance_file_handle_past_comments(fh)
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
    def __str__(self):
        return "\n".join((self.name, self.seq))

class RawEntry(FastQEntry):
    lines = 1
    bowtie_flag = 'r'
    __slots__ = ('name', 'seq')
    def __init__(self, fh):
        self.seq = fh.readline().rstrip('\r\n')

def advance_file_handle_past_comments(fh, out=None):
    pos = fh.tell()
    line = fh.readline()
    while line[0] == "#":
        out and out.write(line)
        pos = fh.tell()
        line = fh.readline()
    fh.seek(pos)


class FastQIndex(object):
    ext = ".fqdx"
    entry_class = FastQEntry

    def is_current(self):
        return is_up_to_date_b(self.filename, self.filename + self.ext)

    @classmethod
    def create(cls, filename, get_next_pos):
        fh = open(filename)
        # iterate past comment lines.
        advance_file_handle_past_comments(fh)

        db = get_db(filename + cls.ext, 'n')
        while True:
            # fh has been moved forward by get_next.
            pos = fh.tell()
            key = get_next_pos(fh)
            if not key: break
            db[key[1:]] = str(pos)
        fh.close()
        db.close()

    def __init__(self, filename):
        self.filename = filename
        self.fh = open(self.filename)
        if not self.is_current():
            try:
                self.create(self.filename, self.entry_class.get_next_position)
            except:
                os.unlink(self.filename + self.ext)
                raise
        self.db = get_db(filename + self.ext, 'r')

    def __getitem__(self, key):
        pos = self.db.get(key, None)
        if pos is None: raise KeyError(key)
        self.fh.seek(long(pos))
        return self.entry_class(self.fh)

    def get_position(self, key):
        return long(self.db.get(key, None))


    def __len__(self):
        return len(self.db)

    def __iter__(self):
        return self.db.iteritems()

    def iterkeys(self):
        return self.db.iterkeys()

    def __contains__(self, key):
        return key in self.db
    contains = __contains__

class FastaIndex(FastQIndex):
    ext = '.fadx'
    entry_class = FastaEntry
    def __init__(self, filename):
        FastQIndex. __init__(self, filename)
