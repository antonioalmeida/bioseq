import unittest
import bioseq as bs

class SequenceTests(unittest.TestCase):

    valid_dna_seq = 'ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA'
    invalid_dna_seq = 'abcd'

    valid_dna = bs.DNASeq(valid_dna_seq)

    def test_constructor(self):
        # Invalid sequence raises exception
        self.assertRaises(AssertionError, bs.DNASeq, self.invalid_dna_seq)

    def test_symbol_frequency(self):
        expected = {'A': 37, 'T': 23, 'G': 17, 'C': 18}
        self.assertDictEqual(self.valid_dna.symbol_frequency(), expected)

    def test_gc_content(self):
        self.assertEqual(self.valid_dna.gc_content(), 0.3684210526315789)

    def test_reverse_complement(self):
        expected = 'TTATTTTATTTATTAACATTCACTATAATGCTGAGTCTGAGTCTGAGCGTAGTCTGATGCGCGATGCTTCAGCTGAGGCTCATTCATAATTTCAT'
        self.assertEqual(self.valid_dna.reverse_complement(), expected)

    def test_transcription(self):
        expected = bs.RNASeq('AUGAAAUUAUGAAUGAGCCUCAGCUGAAGCAUCGCGCAUCAGACUACGCUCAGACUCAGACUCAGCAUUAUAGUGAAUGUUAAUAAAUAAAAUAA')
        self.assertEqual(self.valid_dna.transcription(), expected)

    def test_translation(self):
        expected = bs.AminoacidSeq('MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N')
        self.assertEqual(self.valid_dna.translation(), expected)

    def test_codon_usage(self):
        expected = {'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 1}
        self.assertDictEqual(self.valid_dna.codon_usage('A'), expected)

    def test_reading_frames(self):
        expected = [
            bs.AminoacidSeq('MKL_MSLS_SIAHQTTLRLRLSIIVNVNK_N'),
            bs.AminoacidSeq('_NYE_ASAEASRIRLRSDSDSAL__MLINKI'),
            bs.AminoacidSeq('EIMNEPQLKHRASDYAQTQTQHYSEC__IK_'),
            bs.AminoacidSeq('LFYLLTFTIMLSLSLSVV_CAMLQLRLIHNF'),
            bs.AminoacidSeq('YFIY_HSL_C_V_V_A_SDARCFS_GSFIIS'),
            bs.AminoacidSeq('ILFINIHYNAESESERSLMRDASAEAHS_FH')
        ]

        self.assertListEqual(self.valid_dna.reading_frames(), expected)

    def test_all_proteins_rf(self):
        expected = [
            bs.ProteinSeq('MIEDEMIDMAEIDIAEMD_'),
            bs.ProteinSeq('MIDMAEIDIAEMD_'),
            bs.ProteinSeq('MAEIDIAEMD_'),
            bs.ProteinSeq('MD_'),
        ]
        self.assertListEqual(bs.AminoacidSeq('IAEMIEDEMIDMAEIDIAEMD_OMAED').all_proteins_rf(), expected)

    def test_all_open_reading_frames(self):
        expected = [
            bs.ProteinSeq('MNEPQLKHRASDYAQTQTQHYSEC_'),
            bs.ProteinSeq('MRDASAEAHS_'),
            bs.ProteinSeq('MLSLSLSVV_'),
            bs.ProteinSeq('MSLS_'),
            bs.ProteinSeq('MKL_')
        ]
        self.assertListEqual(self.valid_dna.all_open_reading_frames(), expected)

class AlignmentTests(unittest.TestCase):

    def test_constructor(self):
        alignment = bs.Alignment()

        expected = {'AA': 4, 'AR': -1, 'AN': -2, 'AD': -2, 'AC': 0, 'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -1, 'AL': -1, 'AK': -1, 'AM': -1, 'AF': -2, 'AP': -1, 'AS': 1, 'AT': 0, 'AW': -3, 'AY': -2, 'AV': 0, 'RA': -1, 'RR': 5, 'RN': 0, 'RD': -2, 'RC': -3, 'RQ': 1, 'RE': 0, 'RG': -2, 'RH': 0, 'RI': -3, 'RL': -2, 'RK': 2, 'RM': -1, 'RF': -3, 'RP': -2, 'RS': -1, 'RT': -1, 'RW': -3, 'RY': -2, 'RV': -3, 'NA': -2, 'NR': 0, 'NN': 6, 'ND': 1, 'NC': -3, 'NQ': 0, 'NE': 0, 'NG': 0, 'NH': 1, 'NI': -3, 'NL': -3, 'NK': 0, 'NM': -2, 'NF': -3, 'NP': -2, 'NS': 1, 'NT': 0, 'NW': -4, 'NY': -2, 'NV': -3, 'DA': -2, 'DR': -2, 'DN': 1, 'DD': 6, 'DC': -3, 'DQ': 0, 'DE': 2, 'DG': -1, 'DH': -1, 'DI': -3, 'DL': -4, 'DK': -1, 'DM': -3, 'DF': -3, 'DP': -1, 'DS': 0, 'DT': -1, 'DW': -4, 'DY': -3, 'DV': -3, 'CA': 0, 'CR': -3, 'CN': -3, 'CD': -3, 'CC': 9, 'CQ': -3, 'CE': -4, 'CG': -3, 'CH': -3, 'CI': -1, 'CL': -1, 'CK': -3, 'CM': -1, 'CF': -2, 'CP': -3, 'CS': -1, 'CT': -1, 'CW': -2, 'CY': -2, 'CV': -1, 'QA': -1, 'QR': 1, 'QN': 0, 'QD': 0, 'QC': -3, 'QQ': 5, 'QE': 2, 'QG': -2, 'QH': 0, 'QI': -3, 'QL': -2, 'QK': 1, 'QM': 0, 'QF': -3, 'QP': -1, 'QS': 0, 'QT': -1, 'QW': -2, 'QY': -1, 'QV': -2, 'EA': -1, 'ER': 0, 'EN': 0, 'ED': 2, 'EC': -4, 'EQ': 2, 'EE': 5, 'EG': -2, 'EH': 0, 'EI': -3, 'EL': -3, 'EK': 1, 'EM': -2, 'EF': -3, 'EP': -1, 'ES': 0, 'ET': -1, 'EW': -3, 'EY': -2, 'EV': -2, 'GA': 0, 'GR': -2, 'GN': 0, 'GD': -1, 'GC': -3, 'GQ': -2, 'GE': -2, 'GG': 6, 'GH': -2, 'GI': -4, 'GL': -4, 'GK': -2, 'GM': -3, 'GF': -3, 'GP': -2, 'GS': 0, 'GT': -2, 'GW': -2, 'GY': -3, 'GV': -3, 'HA': -2, 'HR': 0, 'HN': 1, 'HD': -1, 'HC': -3, 'HQ': 0, 'HE': 0, 'HG': -2, 'HH': 8, 'HI': -3, 'HL': -3, 'HK': -1, 'HM': -2, 'HF': -1, 'HP': -2, 'HS': -1, 'HT': -2, 'HW': -2, 'HY': 2, 'HV': -3, 'IA': -1, 'IR': -3, 'IN': -3, 'ID': -3, 'IC': -1, 'IQ': -3, 'IE': -3, 'IG': -4, 'IH': -3, 'II': 4, 'IL': 2, 'IK': -3, 'IM': 1, 'IF': 0, 'IP': -3, 'IS': -2, 'IT': -1, 'IW': -3, 'IY': -1, 'IV': 3, 'LA': -1, 'LR': -2, 'LN': -3, 'LD': -4, 'LC': -1, 'LQ': -2, 'LE': -3, 'LG': -4, 'LH': -3, 'LI': 2, 'LL': 4, 'LK': -2, 'LM': 2, 'LF': 0, 'LP': -3, 'LS': -2, 'LT': -1, 'LW': -2, 'LY': -1, 'LV': 1, 'KA': -1, 'KR': 2, 'KN': 0, 'KD': -1, 'KC': -3, 'KQ': 1, 'KE': 1, 'KG': -2, 'KH': -1, 'KI': -3, 'KL': -2, 'KK': 5, 'KM': -1, 'KF': -3, 'KP': -1, 'KS': 0, 'KT': -1, 'KW': -3, 'KY': -2, 'KV': -2, 'MA': -1, 'MR': -1, 'MN': -2, 'MD': -3, 'MC': -1, 'MQ': 0, 'ME': -2, 'MG': -3, 'MH': -2, 'MI': 1, 'ML': 2, 'MK': -1, 'MM': 5, 'MF': 0, 'MP': -2, 'MS': -1, 'MT': -1, 'MW': -1, 'MY': -1, 'MV': 1, 'FA': -2, 'FR': -3, 'FN': -3, 'FD': -3, 'FC': -2, 'FQ': -3, 'FE': -3, 'FG': -3, 'FH': -1, 'FI': 0, 'FL': 0, 'FK': -3, 'FM': 0, 'FF': 6, 'FP': -4, 'FS': -2, 'FT': -2, 'FW': 1, 'FY': 3, 'FV': -1, 'PA': -1, 'PR': -2, 'PN': -2, 'PD': -1, 'PC': -3, 'PQ': -1, 'PE': -1, 'PG': -2, 'PH': -2, 'PI': -3, 'PL': -3, 'PK': -1, 'PM': -2, 'PF': -4, 'PP': 7, 'PS': -1, 'PT': -1, 'PW': -4, 'PY': -3, 'PV': -2, 'SA': 1, 'SR': -1, 'SN': 1, 'SD': 0, 'SC': -1, 'SQ': 0, 'SE': 0, 'SG': 0, 'SH': -1, 'SI': -2, 'SL': -2, 'SK': 0, 'SM': -1, 'SF': -2, 'SP': -1, 'SS': 4, 'ST': 1, 'SW': -3, 'SY': -2, 'SV': -2, 'TA': 0, 'TR': -1, 'TN': 0, 'TD': -1, 'TC': -1, 'TQ': -1, 'TE': -1, 'TG': -2, 'TH': -2, 'TI': -1, 'TL': -1, 'TK': -1, 'TM': -1, 'TF': -2, 'TP': -1, 'TS': 1, 'TT': 5, 'TW': -2, 'TY': -2, 'TV': 0, 'WA': -3, 'WR': -3, 'WN': -4, 'WD': -4, 'WC': -2, 'WQ': -2, 'WE': -3, 'WG': -2, 'WH': -2, 'WI': -3, 'WL': -2, 'WK': -3, 'WM': -1, 'WF': 1, 'WP': -4, 'WS': -3, 'WT': -2, 'WW': 11, 'WY': 2, 'WV': -3, 'YA': -2, 'YR': -2, 'YN': -2, 'YD': -3, 'YC': -2, 'YQ': -1, 'YE': -2, 'YG': -3, 'YH': 2, 'YI': -1, 'YL': -1, 'YK': -2, 'YM': -1, 'YF': 3, 'YP': -3, 'YS': -2, 'YT': -2, 'YW': 2, 'YY': 7, 'YV': -1, 'VA': 0, 'VR': -3, 'VN': -3, 'VD': -3, 'VC': -1, 'VQ': -2, 'VE': -2, 'VG': -3, 'VH': -3, 'VI': 3, 'VL': 1, 'VK': -2, 'VM': 1, 'VF': -1, 'VP': -2, 'VS': -2, 'VT': 0, 'VW': -3, 'VY': -1, 'VV': 4}
        self.assertDictEqual(alignment.sm, expected)

    def test_recover_global_align(self):
        al = bs.Alignment()
        (s,t) = al.needleman_wunsch('TACT', 'ACTA', -3)

        self.assertListEqual(al.recover_global_align(t, 'TACT', 'ACTA'), ['_TACT', 'ACTA_'])

    def test_global_align_multiple_solutions(self):
        al = bs.Alignment()
        
        (s,t) = al.global_align_multiple_solutions('TACT', 'ACTA', -3)
        self.assertListEqual(s, [[0, -3, -6, -9, -12], [-3, 0, -3, -1, -4], [-6, 1, 0, -3, 3], [-9, -2, 10, 7, 4], [-12, -5, 7, 15, 12]])
        self.assertListEqual(t, [[3, 3, 3, 3, 3], [[2], [1], [3], [1], [3]], [[2], [1], [1], [1, 3], [1]], [[2], [2], [1], [3], [3]], [[2], [2], [2], [1], [3]]])

if __name__ == '__main__':
    unittest.main()

