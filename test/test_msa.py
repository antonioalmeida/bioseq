import unittest
import bioseq as bs

class MSATests(unittest.TestCase):
    sm = bs.alignment.create_substitution_matrix('ATCG', 1, -1)
    msa = bs.msa(["ATAGC", "AACC", "ATGAC"], sm, g=-1) 

    def test(self):
        s1 = "ATAGC"
        s2 = "AACC"
        s3 = "ATGAC"

        sm = bs.alignment.create_substitution_matrix('ATCG', 1, -1)
        (consensus, wrapper) = bs.msa([s1, s2, s3], sm, g=-1).align()
        self.assertEqual(consensus, 'ATGACC')
        self.assertEqual(wrapper.seqs[0], 'AT_AGC')
        self.assertEqual(wrapper.seqs[1], 'A__ACC')
        self.assertEqual(wrapper.seqs[2], 'ATGAC_')

    def test_score_column(self):
        (_, al) = self.msa.align()
        self.assertEqual(self.msa.score_column(al.column(0)), 3)
        self.assertEqual(self.msa.score_column(al.column(1)), -1)
        self.assertEqual(self.msa.score_column(al.column(2)), -2)
        self.assertEqual(self.msa.score_column(al.column(3)), 3)
        self.assertEqual(self.msa.score_column(al.column(4)), -1)
        self.assertEqual(self.msa.score_column(al.column(5)), -1)

    def test_score_sp(self):
        (_, al) = self.msa.align()
        self.assertEqual(self.msa.score_sp(al), 1)
        
    def test_2(self):
        s1 = "PHWAS"
        s2 = "HWASW"
        s3 = "HPHWA"

        msa = bs.msa([s1, s2, s3], g=-8)
        print(msa.align())
