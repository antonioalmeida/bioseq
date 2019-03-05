import unittest
from BioSeq import BioSeq
from DNASeq import DNASeq
from RNASeq import RNASeq
from ProteinSeq import ProteinSeq

class TestDNASeq(unittest.TestCase):

    valid_dna_seq = 'ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA'
    invalid_dna_seq = 'abcd'

    valid_dna = DNASeq(valid_dna_seq)

    def test_constructor(self):
        # Invalid sequence raises exception
        self.assertRaises(AssertionError, DNASeq, self.invalid_dna_seq)

    def test_symbol_frequency(self):
        expected = {'A': 37, 'T': 23, 'G': 17, 'C': 18}
        self.assertDictEqual(self.valid_dna.symbol_frequency(), expected)

    def test_gc_content(self):
        self.assertEqual(self.valid_dna.gc_content(), 0.3684210526315789)

    def test_reverse_complement(self):
        expected = 'TTATTTTATTTATTAACATTCACTATAATGCTGAGTCTGAGTCTGAGCGTAGTCTGATGCGCGATGCTTCAGCTGAGGCTCATTCATAATTTCAT'
        self.assertEqual(self.valid_dna.reverse_complement(), expected)

    def test_transcription(self):
        expected = RNASeq('AUGAAAUUAUGAAUGAGCCUCAGCUGAAGCAUCGCGCAUCAGACUACGCUCAGACUCAGACUCAGCAUUAUAGUGAAUGUUAAUAAAUAAAAUAA')
        self.assertEqual(self.valid_dna.transcription(), expected)

    def test_translation(self):
        expected = ProteinSeq('MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N')
        self.assertEqual(self.valid_dna.translation(), expected)

if __name__ == '__main__':
    unittest.main()
