from bioseq.fasta import read_fasta_file

class Blast:

    db = []
    w = 3
    m = []

    def __init__(self, filename = 'bioseq/res/seqdump.txt', w = 3, fasta=True):
        if fasta:
            self.read_db_fasta(filename)    
        else:
            self.read_db(filename)
        self.w = w

    def read_db_fasta(self, filename):
        self.db = read_fasta_file(filename, 'protein', use_dic=False)

    def read_db(self, filename):
        """From file with sequences line by line read the sequences to a list"""
        with open(filename, 'r') as fd:
            self.db = [line.strip() for line in fd]
        pass

    def db_add_sequence(self, seq):
        """Add an extra sequence to DB"""
        self.db.append(seq)
        
    def build_map(self, query):
        res = {}
        for i in range(len(query)-self.w+1):
            subseq = query[i:i+self.w]
            if subseq in res:
                res[subseq].append(i)
            else:
                res[subseq] = [i]
        
        self.m = res

    def get_hits(self, seq, query):
        res = [] # list of tuples
        for i in range(len(seq)-self.w+1):
            subseq = seq[i:i+self.w]
            if subseq in self.m:
                l = self.m[subseq]
                for ind in l:
                    res.append( (ind,i) )
        return res 

    def extend_hit(self, seq, hit, query):
        stq, sts = hit[0], hit[1]
        matfw = 0       
        k=0
        bestk = 0
        while 2*matfw >= k and stq+self.w+k < len(query) and sts+self.w+k < len(seq):
            if query[stq+self.w+k] == seq[sts+self.w+k]: 
                matfw+=1
                bestk = k+1
            k += 1
        size = self.w + bestk
    
        k = 0
        matbw = 0   
        bestk = 0
        while 2*matbw >= k and stq > k and sts > k:
            if query[stq-k-1] == seq[sts-k-1]: 
                matbw+=1
                bestk = k+1
            k+=1       
        size += bestk
        
        """<index of align. on query, index of align. on sequence, size of align, score>"""
        return (stq-bestk, sts-bestk, size, self.w+matfw+matbw)
        
    def hit_best_score(self, seq, query):
        hits = self.get_hits(seq, query)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extend_hit(seq, h, query)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext

        """<index of align. on query, index of align. on sequence, size of align, score>"""
        return best
 
    def best_alignment(self, q):
        self.build_map(q)
        return self.best_alignment_helper(q)
    
    def best_alignment_helper(self, query, restrictions=[]):
        bestScore = -1.0
        res = (0,0,0,0,0)
        for k in range(0,len(self.db)):
            if k not in restrictions:
                bestSeq = self.hit_best_score(self.db[k], query)
                if bestSeq != ():
                    score = bestSeq[3]  
                    if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                        bestScore = score
                        res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0: return ()
        else: return res
        """<index of align. on query, index of align. on sequence, size of align, score, ???>"""

    def best_n_alignments(self, q, n=10):
        self.build_map(q)
        seqs = []
        species = []

        for k in range(0,len(self.db)):
            if self.db[k].species != q.species and self.db[k].species not in species:
                best = self.hit_best_score(self.db[k], q)
                if best != ():
                    seqs.append(best + (k,))
                    species.append(self.db[k].species)

        scores = sorted(seqs, key=lambda s: s[3], reverse=True)[:n]
        return [ self.db[i[4]] for i in scores], [i[3] for i in scores]





        