import unittest
import bioseq as bs

class MultipleAlignmentTests(unittest.TestCase):

    def test(self):
        s1 = "PHWAS"
        s2 = "HWASW"
        s3 = "HPHWA"

        s1 = "ATAGC"
        s2 = "AACC"
        s3 = "ATGAC"

        sm = bs.alignment.create_substitution_matrix('ATCG', 1, -1)
        msa = bs.msa([s1, s2, s3], sm, g=-1)
        t = msa.align_consensus()
        print(t.seqs)

    def test_2(self):
        s1 = "PHWAS"
        s2 = "HWASW"
        s3 = "HPHWA"

        msa = bs.msa([s1, s2, s3], g=-8)
        t = msa.align_consensus()
        print(t.seqs)
