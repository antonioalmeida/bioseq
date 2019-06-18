import bioseq.alignment as alignment
from bioseq.tree import Tree

class UPGMA:

    def __init__(self, seqs, sm=False):
        self.seqs = seqs
        self.alignment = alignment(sm)
        self.clusters = [str(i) for i,_ in enumerate(seqs)]
        self.trees = {str(i):Tree(seq) for i,seq in enumerate(seqs)}
        self.dist = self._calculate_distances()

    def execute(self):
        t = None

        while(self.dist != [[]]):
            i,j,d = UPGMA._min_distance(self.dist)

            ic= self.clusters[i]; jc= self.clusters[j]
            t1 = self.trees[ic]; t2 = self.trees[jc]

            t = Tree(dist=d/2, left=t1, right=t2)

            UPGMA._update_distances(self.dist, i, j, t1.size()[1], t2.size()[1])

            nc = UPGMA._update_clusters(self.clusters, i, j)
            self.trees[nc] = t          

        return t

    def _update_distances(d, a, b, nA, nB):
            # Swap if the indices are not ordered
            if b < a:
                a, b = b, a

            # For the lower index, reconstruct the entire row (A, i), where i < A
            row = []
            for i in range(0, a):
                dax = d[a][i]*nA; dbx = d[b][i]*nB; nAB = nA+nB
                row.append((dax+dbx)/nAB)
            d[a] = row
            
            # Then, reconstruct the entire column (i, A), where i > A
            #   Note: Since the matrix is lower triangular, row b only contains values for indices < b
            for i in range(a+1, b):
                dax = d[i][a]*nA; dbx = d[b][i]*nB; nAB = nA+nB
                d[i][a] = (dax+dbx)/nAB
                
            #   We get the rest of the values from row i
            for i in range(b+1, len(d)):
                dax = d[i][a]*nA; dbx = d[i][b]*nB; nAB = nA+nB
                d[i][a] = (dax+dbx)/nAB
                # Remove the (now redundant) second index column entry
                del d[i][b]

            # Remove the (now redundant) second index row
            del d[b]

    def _update_clusters(labels, a, b):
        # Swap if the indices are not ordered
        if b < a:
            a, b = b, a

        # Join the labels in the first index
        new_label = "(" + labels[a] + "," + labels[b] + ")"
        labels[a] = new_label

        # Remove the (now redundant) label in the second index
        del labels[b]

        return new_label

    def _calculate_distances(self):
        dist = []
        for j,s1 in enumerate(self.seqs):
            row = []
            for i,s2 in enumerate(self.seqs):
                if i < j:
                    row.append(UPGMA._hamming_distance(s1,s2))
            dist.append(row)
        return dist

    def _hamming_distance(s1, s2):
        d = 0
        for i in range(len(s1)):
            if s1[i] != s2[i]: d += 1
        return d

    def _get_distance(self, i, j):
        for x,y,d in self.distances:
            if (x==i and y==j) or (x==j and y==i): return d
        return d

    def _min_distance(table):
        min_cell = float("inf")
        x, y = -1, -1

        for i in range(len(table)):
            for j in range(len(table[i])):
                if table[i][j] < min_cell:
                    min_cell = table[i][j]
                    x, y = i, j

        # Return the x, y co-ordinate of cell
        return x, y, min_cell


    
