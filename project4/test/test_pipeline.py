import unittest
import bioseq as bs

class PipelineTests(unittest.TestCase):

    fasta_seq = bs.ProteinSeq('MELSVLLFLALLTGLLLLLVQRHPNTHDRLPPGPRPLPLLGNLLQMDRRGLLKSFLRFREKYGDVFTVHLGPRPVVMLCGVEAIREALVDKAEAFSGRGKIAMVDPFFRGYGVIFANGNRWKVLRRFSVTTMRDFGMGKRSVEERIQEEAQCLIEELRKSKGALMDPTFLFQSITANIICSIVFGKRFHYQDQEFLKMLNLFYQTFSLISSVFGQLFELFSGFLKYFPGAHRQVYKNLQEINAYIGHSVEKHRETLDPSAPKDLIDTYLLHMEKEKSNAHSEFSHQNLNLNTLSLFFAGTETTSTTLRYGFLLMLKYPHVAERVYREIEQVIGPHRPPELHDRAKMPYTEAVIYEIQRFSDLLPMGVPHIVTQHTSFRGYIIPKDTEVFLILSTALHDPHYFEKPDAFNPDHFLDANGALKKTEAFIPFSLGKRICLGEGIARAELFLFFTTILQNFSMASPVAPEDIDLTPQECGVGKIPPTYQIRFLPR')

    def test_constructor(self):
        p = bs.pipeline('bioseq/res/source.fasta')
        self.assertEqual(p.seq, self.fasta_seq)

    def test_run_blast(self):
        p = bs.pipeline('bioseq/res/source.fasta', w=10)
        r = p.run_blast()
        print(r)

    def test_run_msa(self):
        p = bs.pipeline('bioseq/res/source.fasta', w=10)
        blast = p.run_blast()
        msa = p.run_msa(blast)

        for seq in msa.seqs:
            print(seq)
