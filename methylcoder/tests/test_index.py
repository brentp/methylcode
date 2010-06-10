import unittest
import os.path as op
from methylcoder.fastqindex import FastQIndex, FastaIndex

PATH = op.dirname(__file__)
DATA = op.join(PATH, "data")

class FastQIndexTest(unittest.TestCase):
    path = op.join(DATA, "sample.fastq")
    def setUp(self):
        pass


    def tearDown(self):
        pass

if __name__ == "__main__":
    unittest.main()
