import unittest
import bioseq.upgma as upgma

class UPGMATests(unittest.TestCase):

    def test_calculate_distances(self):
        s1 = "A_CATATC_AT_"
        s2 = "A_GATATT_AG_"
        s3 = "AACAGATC_T__"
        s4 = "G_CAT__CGATT"
        u = upgma([s1,s2,s3,s4])
        self.assertEqual(u.dist, [[], [3], [4, 6], [5, 8, 9]])        

    def test_execute(self):
        s1 = "A_CATATC_AT_"
        s2 = "A_GATATT_AG_"
        s3 = "AACAGATC_T__"
        s4 = "G_CAT__CGATT"
        u = upgma([s1,s2,s3,s4])
        t = u.execute()
        t.phylo_tree()

    def test_execute_2(self):
        s1 = "AACAGATC_T__"
        s2 = "A_GATATT_AG_"
        s3 = "A_CATATC_AT_"
        s4 = "G_CAT__CGATT"
        u = upgma([s1,s2,s3,s4])
        t = u.execute()
        t.phylo_tree()



