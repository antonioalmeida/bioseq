import bioseq.alignment as alignment
from bioseq.tree import Tree

class UPGMA:

    def __init__(self, seqs, sm=False):
        self.seqs = seqs
        self.alignment = alignment(sm)
        self.i = []
        self.distances = self._calculate_distances()
        self.trees = [Tree(seq) for seq in seqs]

    def execute(self):
        t = None

        while(self.distances):
            (i,j,d) = min(self.distances, key=lambda x:x[2])

            t1 = self.trees[i]; t2 = self.trees[j]
            t = Tree(-1, dist=d/2, left=t1, right=t2)

            self._update_distances(i, j, t1.size()[1], t2.size()[1])
            self.trees[j] = t
        
        return t

    def _update_distances(self, i, j, nA, nB):
        removed = list(filter(lambda x: x[0] != i and x[1] != j, self.distances))

        dist = []
        for ind,(s,x,d) in enumerate(removed):
            if s == j:
                dax = self._get_distance(i,x)
                dbx = self._get_distance(j,x)
                d = (nA*dax + nB*dbx)/(nA+nB)
            dist.append((s,x,d))

        self.distances = dist

    def _calculate_distances(self):
        dist = []
        for i,s1 in enumerate(self.seqs):
            for j,s2 in enumerate(self.seqs):
                if i < j:
                    dist.append((i,j,self._distance(s1,s2)))
        return dist

    def _distance(self, s1, s2):
        d = 0
        for i in range(len(s1)):
            if s1[i] != s2[i]: d += 1
        return d

    def _get_distance(self, i, j):
        for x,y,d in self.distances:
            if x==i and y==j: return d
        return d


    
