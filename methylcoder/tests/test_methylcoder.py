import unittest
import os.path as op
import sys
import os
from methylcoder import bin_paths_from_fasta, _update_conversions
import numpy as np

PATH = op.dirname(__file__)
DATA = op.join(PATH, "data")

class TestBinPaths(unittest.TestCase):
    def setUp(self):
        self.fasta = "/tmp/a.fasta"


    def test_pattern(self):
        paths = bin_paths_from_fasta(self.fasta, pattern_only=True)
        self.assertEquals(paths, 'a.%s.*.bin')

        paths = bin_paths_from_fasta(self.fasta, pattern_only=True, out_dir="/a")
        self.assertEquals(paths, '/a/a.%s.*.bin')

    def test_paths(self):
        paths = bin_paths_from_fasta(self.fasta)
        self.assertEquals(paths, ['a.%s.c.bin', 'a.%s.t.bin', 'a.%s.methyltype.bin'])

class TestCountConversions(unittest.TestCase):
    def setUp(self):
        self.ref = "ACCCTGGA"
        self.aln = "ATTCCGAA"
        self.cc = np.zeros((len(self.ref)), dtype=np.uint32)
        self.tt = np.zeros((len(self.ref)), dtype=np.uint32)

    def test_count(self):
        count = _update_conversions(self.ref, self.aln, 0, "CT", self.cc, self.tt, 10, len(self.ref))
        self.assertEquals(count, 0, (count))

        count = _update_conversions(self.ref, self.aln, 0, "CT", self.cc, self.tt, 0, len(self.ref))
        self.assertEquals(count, 1, (count))

    def test_cts(self):
        _update_conversions(self.ref, self.aln, 0, "CT", self.cc, self.tt, 10, len(self.ref))

        self.assertEquals(self.cc.tolist(),  [0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L])
        self.assertEquals(self.tt.tolist(),  [0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L])

        _update_conversions(self.ref, self.aln, 0, "CT", self.cc, self.tt, 10, len(self.ref))

        self.assertEquals(self.cc.tolist(),  [0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L])
        self.assertEquals(self.tt.tolist(),  [0L, 2L, 2L, 0L, 0L, 0L, 0L, 0L])

    def test_gas(self):
        _update_conversions(self.ref, self.aln, 0, "GA", self.cc, self.tt, 10, len(self.ref))
        for c, t, r, a in zip(self.cc, self.tt, self.ref, self.aln):
            if c != 0:
                self.assertEquals(r, "G")
                self.assertEquals(a, "G")
            if t != 0:
                self.assertEquals(r, "G")
                self.assertEquals(a, "A")


if __name__ == "__main__":
    unittest.main()
