## Graph represented as adjacency list using a dictionary
## keys are vertices
## values of the dictionary represent the list of adjacent vertices of the key node

class Graph:
    """
    Class representing a graph, namely through a dictionary.
    Dictionary keys correspond to vertices.
    Values of the dictionary represent the list of adjacent
    vertices of the correspondent key, i.e. edges.
    """

    def __init__(self, g = {}):
        """Constructor - takes dictionary to fill the graph as input."""
        self.graph = g    

    @staticmethod
    def matrix(m, c=10):
        """Creates a Graph instance from a matrix and a given cut value"""
        g = Graph()

        for j,row in enumerate(m):
            for i,v in enumerate(row):
                if v < c:
                    g.add_edge(i, j) 

        return g

    def print_graph(self):
        """Prints the content of the graph as adjacency list"""
        for v in self.graph.keys():
            print (v, " -> ", self.graph[v])

    def export(self, path='out/network'):
        """Exports graphviz version of the graph"""
        from graphviz import Graph as Vizgraph

        g = Vizgraph('Network', engine='sfdp')
        for i,j in self.get_edges():
            g.edge(str(i), str(j), constraint='false')

        g.render(path, view=False)

    def stats(self):
        """Display information about the graph"""
        nodes = self.get_nodes(); n_nodes = len(nodes)
        edges = self.get_edges(); n_edges = len(edges)
        print('> Graph')
        print('    > %d nodes, %d edges' % (n_nodes, n_edges)); print()

        print('    > Degrees')
        print('        > Highest - ', self.highest_degrees())
        print('        > Mean - %d' % self.mean_degree())
        print('        > Prob - ', list(self.prob_degree().items()))

        print('    > Clustering Coefficients')
        print('        > All - ', list(self.all_clustering_coefs().items()))
        print('        > Highest -', max((self.all_clustering_coefs().items()), key=lambda x:x[1]))
        print('        > Mean - %d' % self.mean_clustering_coef())


    def get_nodes(self):
        """Returns the graph's nodes (keys of dictionary)"""
        return list(self.graph.keys())

    def get_edges(self): 
        """Returns edges in the graph as a list of tuples (origin, destination)"""
        res = []
        for k, v in self.graph.items():
            for dest in v:
                res.append((k, dest))
        return res
      
    def size(self):
        """Returns size of the graph : number of nodes, number of edges"""
        return len(self.graph.keys()), len(self.graph.items())

    def add_vertex(self, v):
        """Add a vertex to the graph; tests if vertex exists not adding if it does"""
        if v not in self.graph:
            self.graph[v] = []

    def add_edge(self, o, d):
        """Add edge to the graph; if vertices do not exist, they are added to the graph"""
        self.add_vertex(o); self.add_vertex(d) # does nothing if vertices exist
        if d not in self.graph[o]:
            self.graph[o].append(d)
        
    def get_successors(self, v):
        """Returns a copy of the successors of node v"""
        return list(self.graph[v])     
             
    def get_predecessors(self, v):
        """Returns the predecessors of node v"""
        res = []
        for k in self.graph.keys(): 
            if v in self.graph[k]: 
                res.append(k)
        return res
    
    def get_adjacents(self, v):
        """Returns the adjacents nodes of node v"""
        suc = self.get_successors(v)
        pred = self.get_predecessors(v)
        res = pred
        for p in suc: 
            if p not in res: res.append(p)
        return res
        
    def out_degree(self, v):
        """Returns the out degree of node v"""
        if v in self.graph:
            return len(self.graph[v])

    def in_degree(self, value):
        """Returns the out degree of node v"""
        count = 0
        for k,v in self.graph.items():
            if k is not value:
               if value in v:
                   count += 1
        return count

    def degree(self, v):
        """Returns the degree (in + out) of node v"""
        return self.out_degree(v) + self.in_degree(v)
        
    def all_degrees(self, deg_type = "inout"):
        """Computes the degree (of a given type) for all nodes. deg_type can be "in", "out", or "inout"""
        assert deg_type in ["inout", "in", "out"], "Invalid degree type."

        degs = {}
        for v in self.graph.keys():
            if deg_type == "out" or deg_type == "inout":
                degs[v] = len(self.graph[v])
            else: degs[v] = 0
        if deg_type == "in" or deg_type == "inout":
            for v in self.graph.keys():
                for d in self.graph[v]:
                    if deg_type == "in" or v not in self.graph[d]:
                        degs[d] = degs[d] + 1
        return degs
    
    def highest_degrees(self, all_deg= None, deg_type = "inout", top= 10):
        """Returns to the top 'top' degrees from the graph"""
        assert deg_type in ["inout", "in", "out"], "Invalid degree type."
        if all_deg is None: 
            all_deg = self.all_degrees(deg_type)

        ord_deg = sorted(list(all_deg.items()), key=lambda x : x[1], reverse = True)
        return list(map(lambda x:x[0], ord_deg[:top]))
        
    def mean_degree(self, deg_type="inout"):
        """Returns the average degree of all nodes: sum of all degrees divided by number of nodes"""
        d = self.all_degrees(deg_type)
        n = len(self.graph.keys())
        return sum(d.values()) / n

    def prob_degree(self, deg_type="inout"):
        """Returns the frequencies of degree values on the graph"""
        degrees = self.all_degrees(deg_type)
        probs = {}
        for _,v in degrees.items():
            if v in probs:
                probs[v] += 1
            else:
                probs[v] = 1

        n = float(len(self.get_nodes()))
        for k in probs.keys():
            probs[k] /= n
        return probs
    
    def clustering_coef(self, v):
        """Returns clustering coefficient for node v"""
        adjs = self.get_adjacents(v)
        if len(adjs) <=1: return 0.0

        # calculate the number of links of the adjacent nodes 
        ligs = 0
        ligs = sum([1 for j in adjs for i in adjs if i != j and i in self.get_adjacents(j)])

        return float(ligs)/(len(adjs)*(len(adjs)-1))
        
    def all_clustering_coefs(self):
        """Returns dictionary with all node's clustering coefficients"""
        return {v: self.clustering_coef(v) for v in self.get_nodes()}

    def mean_clustering_coef(self):
        """Returns the mean value of all node's clustering coefficients"""
        return sum(self.all_clustering_coefs().values()) / self.size()[0] 
