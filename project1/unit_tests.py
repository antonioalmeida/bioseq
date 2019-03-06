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

    def test_codon_usage(self):
        expected = {'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 1}
        self.assertDictEqual(self.valid_dna.codon_usage('A'), expected)

    def test_reading_frames(self):
        expected = [
            ProteinSeq('MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N'),
            ProteinSeq('_NYE_ASAEASRIRLRSDSDSAL__MLINKI'),
            ProteinSeq('EIMNEPQLKHRASDYAQTQTQHYSEC__IK_'),
            ProteinSeq('LFYLLTFTIMLSLSLSVV_CAMLQLRLIHNF'),
            ProteinSeq('YFIY_HSL_C_V_V_A_SDARCFS_GSFIIS'),
            ProteinSeq('ILFINIHYNAESESERSLMRDASAEAHS_FH')
        ]

        self.assertListEqual(self.valid_dna.reading_frames(), expected)



if __name__ == '__main__':
    unittest.main()
