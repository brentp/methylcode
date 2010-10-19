import unittest
from methylcoder import convert_colorspace
from methylcoder.cbowtie import cs2seq, seq2cs

class ColorSpaceTest(unittest.TestCase):

    def test_c2t(self):
        cs = "T11111111111111111111111111111111311"
        out = convert_colorspace(cs, "C", "T")
        self.assertEqual(out, "T11111111111111111111111111111111333")

    def test_cs2seq(self):
        cs = "A321023022"
        self.assertEqual(cs2seq(cs), "ATCAAGCCTC")

    def test_seq2cs(self):
        seq = "ATCAAGCCTC"
        cs = "A321023022"
        self.assertEqual(seq2cs(seq), cs)


if __name__ == "__main__":
    unittest.main()
