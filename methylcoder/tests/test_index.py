import unittest
import os.path as op
import sys
import os
from methylcoder.fastqindex import FastQIndex, FastaIndex, FastQEntry, FastaEntry, guess_index_class
import bsddb

PATH = op.dirname(__file__)
DATA = op.join(PATH, "data")

class GuesserTest(unittest.TestCase):
    def setUp(self):
        self.fastq = op.join(DATA, "sample.fastq.test")
        self.fasta = op.join(DATA, "sample.fasta.test")

    def test_fastq(self):
        self.assertEquals(FastQIndex, guess_index_class(self.fastq))

    def test_fasta(self):
        self.assertEquals(FastaIndex, guess_index_class(self.fasta))

    def test_bad(self):
        self.assertRaises(AssertionError, guess_index_class, self.fasta + "ASDF")


class FastQIndexTest(unittest.TestCase):
    base_file = "sample.fastq.test"
    klass = FastQIndex
    header_start = "@"
    def setUp(self):
        self.path = op.join(DATA, self.base_file)
        self.idx_path = self.path + self.klass.ext

    def test_create(self):
        self.assert_(op.exists(self.path))
        fi = self.klass(self.path)
        self.assert_(op.exists(self.idx_path))

    def test_len(self):
        fi = self.klass(self.path)
        nlines = sum(1 for i in open(self.path))
        self.assertEqual(nlines / self.klass.entry_class.lines, len(fi))

    def test_contains(self):
        fi = self.klass(self.path)
        for header in (line.strip() for line in open(self.path) \
                                       if line[0] == self.header_start):

            self.assert_(header[1:] in fi, (header, fi.iterkeys().next()))

    def test_sequence(self):
        fi = self.klass(self.path)
        key, pos = iter(fi).next()
        obj = fi[key]
        pos = int(pos)
        fh = open(self.path, "r")
        fh.seek(pos)
        entry = self.klass.entry_class(fh)
        self.assertEquals(obj.seq, entry.seq, (obj, entry, pos))

    def tearDown(self):
        os.unlink(self.idx_path)


class FastaIndexTest(FastQIndexTest):
    base_file = "sample.fasta.test"
    klass = FastaIndex
    header_start = ">"

if __name__ == "__main__":
    unittest.main()
