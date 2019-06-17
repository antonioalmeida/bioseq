import bioseq.alignment as alignment
from collections import Counter

class MultipleAlignment():

    def __init__(self, seqs, sm=False, g=-1):
        self.seqs = seqs # list of BioSeq objects
        self.alignment = alignment(sm) # an Alignment instance
        self.g = g
    
    def n_seqs(self):
        return len(self.seqs)
    
    def add_seq(self, al, seq):
        res = []
        for i in range(len(al.seqs)+1): res.append('')
        cons = al.consensus()
        print('CONSENSUS')
        print(cons)

        new_al = self.alignment.run_global_align_multiple_solutions(cons, seq, g=self.g, debug=True)[0]
        print('NEW AL')
        print(new_al)
        new_al = AlignmentWrapper(new_al)

        orig = 0
        for i in range(len(new_al)):
            print('i - ' + str(i))
            if new_al[0,i] == '_':
                for k in range(len(al.seqs)):
                    res[k] += "_"
            else:
                for k in range(len(al.seqs)):
                    print('k - ' + str(k))
                    res[k] += al[k, orig]
                orig += 1
        res[len(al.seqs)] = new_al.seqs[1]

        return(AlignmentWrapper(res))

    def align_consensus(self):
        [seq1, seq2] = self.seqs[0:2]
        al = self.alignment.run_global_align_multiple_solutions(seq1, seq2, g=self.g, debug=True)[0]

        al = AlignmentWrapper(al)

        res = al
        for i in range(2, self.n_seqs()):
            res = self.add_seq(res, self.seqs[i])

        return res
        
    def ScoreColumn(self, charsCol):
        pass
        
    def scoreSP (self, alignment):
        pass
   
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
