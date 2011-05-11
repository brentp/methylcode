import unittest
import os.path as op
import sys
import os
from methylcoder import bin_paths_from_fasta, _update_conversions, methylation_summary
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
        count = _update_conversions(self.ref, self.aln, 0, "CT", self.cc, self.tt, 10, len(self.ref), "")
        self.assertEquals(count, 0, (count))

        count = _update_conversions(self.ref, self.aln, 0, "CT", self.cc, self.tt, 0, len(self.ref), "")
        self.assertEquals(count, 1, (count))

    def test_cts(self):
        _update_conversions(self.ref, self.aln, 0, "CT", self.cc, self.tt, 10, len(self.ref), "")

        self.assertEquals(self.cc.tolist(),  [0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L])
        self.assertEquals(self.tt.tolist(),  [0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L])

        _update_conversions(self.ref, self.aln, 0, "CT", self.cc, self.tt, 10, len(self.ref), "")

        self.assertEquals(self.cc.tolist(),  [0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L])
        self.assertEquals(self.tt.tolist(),  [0L, 2L, 2L, 0L, 0L, 0L, 0L, 0L])

    def test_gas(self):
        _update_conversions(self.ref, self.aln, 0, "GA", self.cc, self.tt, 10, len(self.ref), "")
        for c, t, r, a in zip(self.cc, self.tt, self.ref, self.aln):
            if c != 0:
                self.assertEquals(r, "G")
                self.assertEquals(a, "G")
            if t != 0:
                self.assertEquals(r, "G")
                self.assertEquals(a, "A")

class TestMethylationSummary(unittest.TestCase):
    def setUp(self):
        self.cs = np.zeros((10, ), 'i')
        self.ts = np.zeros((10, ), 'i')
        self.cs[4:6] = 2
        self.ts[4:6] = 2

    def test_cg_ctx(self):
        CTX = 'CG'
        mt = np.zeros((10, ), 'i')

        # only set CG
        mt[4:6] = 1
        m = methylation_summary(self.cs, self.ts, mt)
        tsum = self.ts.sum()
        csum = self.cs.sum()
        self.assertEquals(m[CTX]['cs'], csum)
        self.assertEquals(m[CTX]['ts'], tsum)
        self.assertEquals(m[CTX]['methylation'], csum / float(csum + tsum))

    def test_chg_ctx(self):
        CTX = 'CHG'
        mt = np.zeros((10, ), 'i')

        # only set CG
        mt[4:6] = 2
        mt[1:5] = 5
        m = methylation_summary(self.cs, self.ts, mt)
        tsum = self.ts.sum()
        csum = self.cs.sum()
        self.assertEquals(m[CTX]['cs'], csum)
        self.assertEquals(m[CTX]['ts'], tsum)
        self.assertEquals(m[CTX]['methylation'], csum / float(csum + tsum))

    def test_chh_ctx(self):
        CTX = 'CHH'
        mt = np.zeros((10, ), 'i')

        # only set CG
        mt[4:6] = 3
        mt[1:5] = 6
        m = methylation_summary(self.cs, self.ts, mt)
        tsum = self.ts.sum()
        csum = self.cs.sum()
        self.assertEquals(m[CTX]['cs'], csum)
        self.assertEquals(m[CTX]['ts'], tsum)
        self.assertEquals(m[CTX]['methylation'], csum / float(csum + tsum))


if __name__ == "__main__":
    unittest.main()
