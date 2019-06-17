import unittest
import bioseq.upgma as upgma

class UPGMATests(unittest.TestCase):

    def test_calculate_distances(self):
        s1 = "A_CATATC_AT_"
        s2 = "A_GATATT_AG_"
        s3 = "AACAGATC_T__"
        s4 = "G_CAT__CGATT"
        u = upgma([s1,s2,s3,s4])
        self.assertEqual(u.distances, [(0, 1, 3), (0, 2, 4), (0, 3, 5), (1, 2, 6), (1, 3, 8), (2, 3, 9)])        

    def test_execute(self):
        s1 = "A_CATATC_AT_"
        s2 = "A_GATATT_AG_"
        s3 = "AACAGATC_T__"
        s4 = "G_CAT__CGATT"
        u = upgma([s1,s2,s3,s4])
        t = u.execute()
        t.print_tree()

