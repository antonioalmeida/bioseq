from BioSeq import BioSeq
from SequenceType import SequenceType

class DNASeq(BioSeq):

    nucleotides_dic = set(['A', 'T', 'C', 'G'])
    
    def __init__(self, seq):
        super().__init__(seq, SequenceType.DNA)

    def is_valid(self):
        return all(i in DNASeq.nucleotides_dic for i in self.seq)

    def gc_content(self):
        symbols = ['G', 'C']
        count = 0

        for i in self.seq:
            if i in symbols:
                count += 1
        return count/len(self.seq)
