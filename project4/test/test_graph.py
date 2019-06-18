import unittest
import bioseq as bs

class GraphTests(unittest.TestCase):

    def test_vertices(self):
        g = bs.graph.Graph()
        g.add_vertex(1)
        g.add_vertex(2)
        g.add_vertex(3)
        g.add_vertex(4)
        self.assertEqual(len(g.get_nodes()), 4)

    def test_edges(self):
        g = bs.graph.Graph()
        g.add_vertex(1)
        g.add_vertex(2)
        g.add_vertex(3)
        g.add_vertex(4)

        g.add_edge(1, 2)
        g.add_edge(2, 3)
        g.add_edge(3, 2)
        g.add_edge(3, 4)
        g.add_edge(4, 2)

        self.assertEqual(len(g.get_edges()), 5)

    def test_degrees(self):
        g = bs.graph.Graph()
        g.add_vertex(1)
        g.add_vertex(2)
        g.add_vertex(3)
        g.add_vertex(4)

        g.add_edge(1, 2)
        g.add_edge(2, 3)
        g.add_edge(3, 2)
        g.add_edge(3, 4)
        g.add_edge(4, 2)

        self.assertEqual(g.in_degree(2), 3)
        self.assertEqual(g.out_degree(2), 1)
        self.assertEqual(g.degree(2), 4)

    def test_matrix(self):
        m = [ [],                         
              [19],                       
              [27, 31],                   
              [8, 18, 26],                
              [33, 36, 41, 31],           
              [18, 1, 32, 17, 35],        
              [13, 13, 29, 14, 28, 12]    
            ]

        g2 = bs.graph.Graph(m=m, c=15)
        self.assertEqual(g2.get_nodes(), [0, 3, 1, 5, 6])
        self.assertEqual(g2.get_edges(), [(0, 3), (0, 6), (3, 6), (1, 5), (1, 6), (5, 6)])


  
