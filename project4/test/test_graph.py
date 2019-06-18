import unittest
import bioseq as bs

class GraphTests(unittest.TestCase):

    def test_matrix(self):
        m = [ [],                         
              [19],                       
              [27, 31],                   
              [8, 18, 26],                
              [33, 36, 41, 31],           
              [18, 1, 32, 17, 35],        
              [13, 13, 29, 14, 28, 12]    
            ]

        g = bs.Graph.matrix(m, 15)
        self.assertEqual(g.get_nodes(), [0, 3, 1, 5, 6])
        self.assertEqual(g.get_edges(), [(0, 3), (0, 6), (3, 6), (1, 5), (1, 6), (5, 6)])

        g.stats()

    def test_stats(self):
        g = bs.Graph()
        g.export()
