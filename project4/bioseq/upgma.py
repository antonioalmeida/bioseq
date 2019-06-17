import bioseq.alignment as alignment

class UPGMA:

    def __init__(self, seqs, sm=False):
        self.seqs = seqs
        self.alignment = alignment(sm)
        self.distances = self._calculate_distances()

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

    
