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
        self.blast = p.run_blast()

    def test_run_msa(self):
        p = bs.pipeline('bioseq/res/source.fasta')
        b = p.run_blast()
        p.run_msa(b)

