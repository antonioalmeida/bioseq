import unittest
import bioseq as bs

class BlastTests(unittest.TestCase):

    expected_db = ['atatagtgctagatctgatcgatgctaaatgctagatagtgatcgtagagctgatagcgtagatatagatgcgctgctagaatgctgatagctgaaatagatcgatagtcgatagtcgctgatagtctgtagctgatgcgctgatggatgtttgggtacacgatgctgatacgcatgatgatgatgaaca', 'atatagtctaatctgatcgcatgctaatgcgtagatagtgatcgtagagtgatagcgtaggatataatgcgctgctagaatctgataagctgaaatagatcgatagtcatagtcgctgatagtctgtaagtgatgcgcgatggatgtttggtactacgatgctatagctcgagatcgaagcgatgct', 'ataaatagtgctagatctgatcactcgtcgtcgctagatagtgattcgtagagctgatagcgtagatatagatgcgctgctagaatgctgatagctgaaatactcgcgatgctgctgatagtctgtagctgatgcgctgatggatgtttgggtacacgatgctgactgcgatatgagatcgatagcg', 'atagctcgatgcttagatctcgcgtatgctgctagataagagctgctgagctgatcggatgcctcgcgctcgcgcgctgaggctcggatagctagctgagcgctcgatagcgcgttcgctggatcgcgtatagcgctgaagctcccggctagctgtctgtaaatcggatctgatctcgctctatact', 'cgacgacgacgacgaatgatgatgatgaccgccgtagagctagctgagctgcttagctgatcgcgatatgagagagagctatagagcttttttttttttgggggggggggggaaaaaaaaaaaggggggggatagagatagctgaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa']

    def test_constructor(self):
        b = bs.blast('bioseq/res/seqBlast.txt', fasta=False)
        self.assertEqual(b.w, 3)
        self.assertEqual(b.db, self.expected_db)

    def test_fasta_constructor(self):
        b = bs.blast()
        self.assertEqual(b.w, 3)
        self.assertEqual(len(b.db), 21)

    def test_1(self):
        b = bs.blast('bioseq/res/seqBlast.txt', w=11, fasta=False)
        q = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
        r = b.best_alignment(q)
        self.assertEqual(r, (1, 38, 149, 108, 3))
        
    def test_2(self):
        b = bs.blast('bioseq/res/seqBlast.txt', w=11, fasta=False)
        q = "cgacgacgacgacgaatgatg"
        r = b.best_alignment(q)
        self.assertEqual(r, (0, 0, 21, 21, 4))

    def test_seqdump(self):
        s = 'MELSVLLFLALLTGLLLLLVQRHPNTHDRLPPGPRPLPLLGNLLQMDRRGLLKSFLRFREKYGDVFTVHLGPRPVVMLCGVEAIREALVDKAEAFSGRGKIAMVDPFFRGYGVIFANGNRWKVLRRFSVTTMRDFGMGKRSVEERIQEEAQCLIEELRKSKGALMDPTFLFQSITANIICSIVFGKRFHYQDQEFLKMLNLFYQTFSLISSVFGQLFELFSGFLKYFPGAHRQVYKNLQEINAYIGHSVEKHRETLDPSAPKDLIDTYLLHMEKEKSNAHSEFSHQNLNLNTLSLFFAGTETTSTTLRYGFLLMLKYPHVAERVYREIEQVIGPHRPPELHDRAKMPYTEAVIYEIQRFSDLLPMGVPHIVTQHTSFRGYIIPKDTEVFLILSTALHDPHYFEKPDAFNPDHFLDANGALKKTEAFIPFSLGKRICLGEGIARAELFLFFTTILQNFSMASPVAPEDIDLTPQECGVGKIPPTYQIRFLPR'
        seq = bs.ProteinSeq(s, 'NP_000758.1 cytochrome P450 2B6 precursor [Homo sapiens]')

        b = bs.blast('bioseq/res/seqdump.txt', 3)
        seqs,scores = b.best_n_alignments(seq)
        self.assertEqual(scores, [484, 484, 479, 477, 453, 451, 451, 450, 449, 449])

        

