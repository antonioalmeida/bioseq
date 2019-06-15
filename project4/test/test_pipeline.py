import unittest
import bioseq as bs

class PipelineTests(unittest.TestCase):

    def test_constructor(self):
        test = bs.pipeline('bioseq/res/source.fasta', 'protein')
        print(test.seq)