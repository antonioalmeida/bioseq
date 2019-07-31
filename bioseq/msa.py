import bioseq.alignment as alignment
from bioseq.proteinseq import ProteinSeq
from collections import Counter
from itertools import combinations

class MultipleAlignment():

    def __init__(self, seqs, sm=False, g=-3, species=False):
        self.seqs = seqs # list of BioSeq objects
        self.alignment = alignment(sm) # an Alignment instance
        self.g = g
        self.species = []

        if species:
            self.species = [s.species for s in self.seqs]
    
    def n_seqs(self):
        return len(self.seqs)
    
    def add_seq(self, al, seq):
        cons = al.consensus()
        new_al = self.alignment.run_global_align_multiple_solutions(cons, seq, g=self.g)[0]
        new_al = AlignmentWrapper(new_al)

        return self._normalize(al, new_al)

    def align(self):
        [seq1, seq2] = self.seqs[0:2]
        al = self.alignment.run_global_align_multiple_solutions(seq1, seq2, g=self.g)[0]

        al = AlignmentWrapper(al)

        res = al
        for i in range(2, self.n_seqs()):
            res = self.add_seq(res, self.seqs[i])

        final = []
        if self.species:
            for i,seq in enumerate(res):
                final.append(ProteinSeq(seq, self.species[i]))
            final = AlignmentWrapper(final)
        else:
            final = res
                
        return (final.consensus(), final)
        
    def score_column(self, col):
        score = 0

        filter(lambda x: x != '_', col)
        pairs = list(combinations(col, 2))

        for x,y in pairs:
            score += self._score_pair(x,y)
        return score
        
    def score_sp(self, alignment):
        score = 0
        for i in range(len(alignment)):
            score += self.score_column(alignment.column(i))

        return score

    def _score_pair(self, c1, c2):
        # source for definition - https://www.cs.princeton.edu/~mona/Lecture/msa1.pdf
        if c1 == '_' and c2 == '_': return 0
        elif c1 == '_' or c2 == '_': return self.g
        else: return self.alignment.sm[c1+c2]

    def _normalize(self, previous, new):
        res = []; n_seqs = len(previous.seqs)
        for _ in range(n_seqs+1): res.append('')

        o = 0
        for i in range(len(new)):
            if new[0,i] == '_':
                for k in range(n_seqs):
                    res[k] += "_"
            else:
                for k in range(n_seqs):
                    res[k] += previous[k, o]
                o += 1

        res[n_seqs] = new.seqs[1]
        return AlignmentWrapper(res)

class AlignmentWrapper():

    def __init__(self, seqs, seq_type = "protein"):
        self.seqs = seqs
        self.seq_type = seq_type

    def __len__(self): # number of columns
        return len(self.seqs[0])
    
    def __getitem__(self, n):
        if type(n) is tuple and len(n) ==2: 
            i, j = n
            return self.seqs[i][j]
        elif type(n) is int: return self.seqs[n]
        return None
    
    def __str__(self):
        res = ""
        for seq in self.seqs:
            res += "\n" + seq 
        return res
    
    def num_seqs(self):
        return len(self.seqs)
       
    def column (self, indice):
        res = []
        for k in range(len(self.seqs)):
            res.append(self.seqs[k][indice])
        return res
    
    def consensus(self):
        res = ''
        for i in range(len(self)):
            col = [self.seqs[j][i] for j in range(len(self.seqs))]

            dic = {}
            for j in col:
                if j in dic:
                    dic[j] += 1
                elif j != '_':
                    dic[j] = 1
            
            max_value = max(dic)
            c = min(dic.keys(), key=lambda x: x == max_value)
            res += c

        return res

    def print(self):
        print('> Visualizing MSA')

        for seq in self.seqs:
            if type(seq) is str:
                print('> ' + seq)
            else:
                print('> ' + seq.seq + str(seq))

        print()
        print('    > Consensus')
        print('    > ' + self.consensus())
