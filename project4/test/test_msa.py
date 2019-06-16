import unittest
import bioseq as bs

class MultipleAlignmentTests(unittest.TestCase):

    def test(self):
        s1 = "PHWAS"
        s2 = "HWASW"
        s3 = "HPHWA"
        msa = bs.msa([s1, s2, s3])
        t = msa.align_consensus()
        print(t.seqs)
