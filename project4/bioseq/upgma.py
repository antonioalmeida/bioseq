import bioseq.alignment as alignment
from bioseq.tree import Tree

class UPGMA:

    def __init__(self, seqs, sm=False):
        self.seqs = seqs
        self.alignment = alignment(sm)
        self.i = []
        self.distances = self._calculate_distances()

        self.labels = [str(i) for i,_ in enumerate(seqs)]
        #self.trees = [Tree(seq) for seq in seqs]

        self.trees = {str(i):Tree(seq) for i,seq in enumerate(seqs)}

        self.dist = self._calculate_distances_2()

    def execute(self):
        t = None

        print(self.dist)
        while(self.dist != [[]]):
            #(i,j,d) = min(self.distances, key=lambda x:x[2])
            #print('idj ' + str(i) +' - '+str(j)+' - '+str(d))

            i,j,d = UPGMA._min_distance(self.dist)

            i_cluster = self.labels[i];
            j_cluster = self.labels[j];
            t1 = self.trees[i_cluster]; t2 = self.trees[j_cluster]

            t = Tree(dist=d/2, left=t1, right=t2)
            #self._update_distances(i, j, t1.size()[1], t2.size()[1])

            print('IDJ - ', i,j,d)
            UPGMA._join_table(self.dist, i, j, t1.size()[1], t2.size()[1])

            new_cluster = UPGMA._join_labels(self.labels, i, j)
 
            #del self.trees[i]
            #del self.trees[j]

            self.trees[new_cluster] = t          
            print(self.dist)

        for t,y in self.trees.items():
            y.print_tree_2()
        
        return t

    def _join_table(d, a, b, nA, nB):
            # Swap if the indices are not ordered
            if b < a:
                a, b = b, a

            # For the lower index, reconstruct the entire row (A, i), where i < A
            row = []
            for i in range(0, a):
                dax = d[a][i]*nA; dbx = d[b][i]*nB; nAB = nA+nB

                print('DAX - DBX ', (dax, dbx))
                row.append((dax+dbx)/nAB)
            d[a] = row
            
            # Then, reconstruct the entire column (i, A), where i > A
            #   Note: Since the matrix is lower triangular, row b only contains values for indices < b
            for i in range(a+1, b):
                dax = d[i][a]*nA; dbx = d[b][i]*nB; nAB = nA+nB
                print('DAX - DBX ', (dax, dbx))
                d[i][a] = (dax+dbx)/nAB
                
            #   We get the rest of the values from row i
            for i in range(b+1, len(d)):
                print('dia -',  d[i][a], nA)
                dax = d[i][a]*nA; dbx = d[i][b]*nB; nAB = nA+nB
                print('DAX - DBX ', (dax, dbx))
                print('D')
                print(nAB)
                d[i][a] = (dax+dbx)/nAB
                # Remove the (now redundant) second index column entry
                del d[i][b]

            # Remove the (now redundant) second index row
            del d[b]

    def _join_labels(labels, a, b):
        # Swap if the indices are not ordered
        if b < a:
            a, b = b, a

        # Join the labels in the first index
        new_label = "(" + labels[a] + "," + labels[b] + ")"
        labels[a] = new_label

        # Remove the (now redundant) label in the second index
        del labels[b]

        return new_label


    def _update_distances(self, i, j, nA, nB):
        removed = list(filter(lambda x: x[0] == i or x[1] == i, self.distances))

        dist = []
        new_i = len(self.trees)-1
        for ind,(x,y,d) in enumerate(removed):
            if x == j:  
                print('TEST')
                dax = self._get_distance(i,y)
                dbx = self._get_distance(j,y)

                print('i,j - %i - %i', (i, j))
                print('dax - %i - %i', (dax, dbx))
                d = (nA*dax + nB*dbx)/(nA+nB) 
                dist.append((new_i,y,d))
            elif j == i:
                dax = self._get_distance(i,x)
                dbx = self._get_distance(j,x)

                print('i,j - %i - %i', (i, j))
                print('dax - %i - %i', (dax, dbx))
                d = (nA*dax + nB*dbx)/(nA+nB) 
                dist.append((x,new_i,d))
            else:
                dist.append((x,y,d))

        self.distances = dist

    def _calculate_distances_2(self):
        dist = []
        for j,s1 in enumerate(self.seqs):
            row = []
            for i,s2 in enumerate(self.seqs):
                if i < j:
                    row.append(UPGMA._hamming_distance(s1,s2))
            dist.append(row)
        return dist

    def _calculate_distances(self):
        dist = []
        for i,s1 in enumerate(self.seqs):
            for j,s2 in enumerate(self.seqs):
                if i < j:
                    dist.append((i,j,UPGMA._hamming_distance(s1,s2)))
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


    
