import unittest
import bioseq as bs

class BlastTests(unittest.TestCase):

    expected_db = ['atatagtgctagatctgatcgatgctaaatgctagatagtgatcgtagagctgatagcgtagatatagatgcgctgctagaatgctgatagctgaaatagatcgatagtcgatagtcgctgatagtctgtagctgatgcgctgatggatgtttgggtacacgatgctgatacgcatgatgatgatgaaca', 'atatagtctaatctgatcgcatgctaatgcgtagatagtgatcgtagagtgatagcgtaggatataatgcgctgctagaatctgataagctgaaatagatcgatagtcatagtcgctgatagtctgtaagtgatgcgcgatggatgtttggtactacgatgctatagctcgagatcgaagcgatgct', 'ataaatagtgctagatctgatcactcgtcgtcgctagatagtgattcgtagagctgatagcgtagatatagatgcgctgctagaatgctgatagctgaaatactcgcgatgctgctgatagtctgtagctgatgcgctgatggatgtttgggtacacgatgctgactgcgatatgagatcgatagcg', 'atagctcgatgcttagatctcgcgtatgctgctagataagagctgctgagctgatcggatgcctcgcgctcgcgcgctgaggctcggatagctagctgagcgctcgatagcgcgttcgctggatcgcgtatagcgctgaagctcccggctagctgtctgtaaatcggatctgatctcgctctatact', 'cgacgacgacgacgaatgatgatgatgaccgccgtagagctagctgagctgcttagctgatcgcgatatgagagagagctatagagcttttttttttttgggggggggggggaaaaaaaaaaaggggggggatagagatagctgaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa']

    def test_constructor(self):
        b = bs.blast('bioseq/res/seqBlast.txt')
        self.assertEqual(b.w, 3)
        self.assertEqual(b.db, self.expected_db)

    def test_1(self):
        b = bs.blast('bioseq/res/seqBlast.txt', 11)
        q = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
        r = b.bestAlignment(q)
        self.assertEqual(r, (1, 38, 149, 108, 3))
        
    def test_2(self):
        b = bs.blast('bioseq/res/seqBlast.txt', 11)
        q = "cgacgacgacgacgaatgatg"
        r = b.bestAlignment(q)
        self.assertEqual(r, (0, 0, 21, 21, 4))
        

