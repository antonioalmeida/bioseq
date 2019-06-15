# -*- coding: utf-8 -*-

class Blast:

    db = []
    w = 3
    m = []

    def __init__(self, filename = 'bioseq/res/seqdump.txt', w = 3):
        self.read_db(filename)
        self.w = w

    def read_db(self, filename):
        """From file with sequences line by line read the sequences to a list"""
        with open(filename, 'r') as fd:
            self.db = [line.strip() for line in fd]
        pass

    def addSequenceDB(self, seq):
        """Add an extra sequence to DB"""
        self.db.append(seq)
        
    def buildMap (self, query):
        res = {}
        for i in range(len(query)-self.w+1):
            subseq = query[i:i+self.w]
            if subseq in res:
                res[subseq].append(i)
            else:
                res[subseq] = [i]
        
        self.m = res

    def getHits(self, seq, query):
        res = [] # list of tuples
        for i in range(len(seq)-self.w+1):
            subseq = seq[i:i+self.w]
            if subseq in self.m:
                l = self.m[subseq]
                for ind in l:
                    res.append( (ind,i) )
        return res 

    def extendsHit (self, seq, hit, query):
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
        
        return (stq-bestk, sts-bestk, size, self.w+matfw+matbw)
        
    def hitBestScore(self, seq, query):
        hits = self.getHits(seq, query)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extendsHit(seq, h, query)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext
        return best
 
 
    def bestAlignment (self, query):
        self.buildMap(query)
        bestScore = -1.0
        res = (0,0,0,0,0)
        for k in range(0,len(self.db)):
            bestSeq = self.hitBestScore(self.db[k], query)
            if bestSeq != ():
                score = bestSeq[3]  
                if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                    bestScore = score
                    res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0: return ()
        else: return res

def test1():
    mb = Blast("bioseq/res/seqBlast.txt", 11)
    print(mb.db)
    query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
    r = mb.bestAlignment(query)
    print(r)


def test2():
    mb = Blast("bioseq/res/seqBlast.txt", 11)
    query2 = "cgacgacgacgacgaatgatg"
    r = mb.bestAlignment(query2)
    print(r)       

if __name__ == '__main__':
    test1()
    test2()
