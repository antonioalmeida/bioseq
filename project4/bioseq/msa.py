import bioseq.alignment 

class MultipleAlignment():

    def __init__(self, seqs, alignment=bioseq.alignment()):
        self.seqs = seqs # list of BioSeq objects
        self.alignment = alignment # an Alignment instance
    
    def num_seqs(self):
        return len(self.seqs)
    
    def add_seq_alignment (self, alignment, seq):
        pass
        # res = []
        # for i in range(len(alignment.seqs)+1):
        #     res.append("")
        # # create consensus from give alignments
        # cons = BioSeq(alignment.consensus(), alignment.seq_type)

        # for i in range(len(align2)):
        #     if align2[0,i]== '-':
        #         for k in range(len(alignment.seqs)):
        #             res[k] += "-"
        #     else:
        #         for k in range(len(alignment.seqs)):
        #             res[k] += alignment[k,orig]
        #         orig+=1
        # res[len(alignment.seqs)] = align2.seqs[1]

        # return AlignmentWrapper(res, alignment.seq_type) 
    
    def align_consensus(self):
        [seq1, seq2] = self.seqs[0:2]
        initial_align = self.alignment.run_global_align_multiple_solutions(seq1, seq2, g=-8, debug=True)
        print(initial_align)
        
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
    
    def consensus (self):
        cons = ""
        for i in range(len(self)):
            cont = {}
            for k in range(len(self.seqs)):
                c = self.seqs[k][i]
                if c in cont:
                    cont[c] = cont[c] + 1
                else:
                    cont[c] = 1
                maximum = 0
                cmax = None
            for ke in cont.keys():
                if ke != "-" and cont[ke] > maximum:
                    maximum = cont[ke]
                    cmax = ke
            cons = cons + cmax
        return cons 


if __name__ == "__main__":   
    s1 = "PHWAS"
    s2 = "HWASW"
    s3 = "HPHWA"
    msa = MultipleAlignment([s1, s2, s3])
    msa.align_consensus()