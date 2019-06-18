import unittest
import bioseq as bs

class PipelineTests(unittest.TestCase):

    fasta_seq = bs.ProteinSeq('MELSVLLFLALLTGLLLLLVQRHPNTHDRLPPGPRPLPLLGNLLQMDRRGLLKSFLRFREKYGDVFTVHLGPRPVVMLCGVEAIREALVDKAEAFSGRGKIAMVDPFFRGYGVIFANGNRWKVLRRFSVTTMRDFGMGKRSVEERIQEEAQCLIEELRKSKGALMDPTFLFQSITANIICSIVFGKRFHYQDQEFLKMLNLFYQTFSLISSVFGQLFELFSGFLKYFPGAHRQVYKNLQEINAYIGHSVEKHRETLDPSAPKDLIDTYLLHMEKEKSNAHSEFSHQNLNLNTLSLFFAGTETTSTTLRYGFLLMLKYPHVAERVYREIEQVIGPHRPPELHDRAKMPYTEAVIYEIQRFSDLLPMGVPHIVTQHTSFRGYIIPKDTEVFLILSTALHDPHYFEKPDAFNPDHFLDANGALKKTEAFIPFSLGKRICLGEGIARAELFLFFTTILQNFSMASPVAPEDIDLTPQECGVGKIPPTYQIRFLPR')

    pipeline = None
    blast = None

    def test_constructor(self):
        p = bs.pipeline('bioseq/res/source.fasta')
        self.assertEqual(p.seq, self.fasta_seq)

    def test_run_blast(self):
        p = bs.pipeline('bioseq/res/source.fasta')
        p.run_blast()
        p.run_msa()
        p.run_upgma()

        p.visualize_msa()
        p.phylo()
        p.text()

    # def test_run_msa(self):
    #     p = bs.pipeline('bioseq/res/source.fasta')
    #     sm = bs.alignment.create_substitution_matrix('ABCDEFGHIJKLMNOPQRSTUVWXYZ', 1, -1)
    #     b = p.run_blast()
    #     p.run_msa(b, sm)

    # def test_run_upga(self):
    #     p = bs.pipeline('bioseq/res/source.fasta')
    #     sm = bs.alignment.create_substitution_matrix('ABCDEFGHIJKLMNOPQRSTUVWXYZ', 1, -1)
    #     b = p.run_blast()
    #     al = p.run_msa(b)
    #     p.run_upgma(al)
    #     p.phylo()
    #     p.text()



