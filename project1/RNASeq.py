from BioSeq import BioSeq
from SequenceType import SequenceType

class RNASeq(BioSeq):

    nucleotides_dic = set(['A', 'U', 'C', 'G'])

    complement = {"A": "U",
                  "U": "A",
                  "G": "C",
                  "C": "G"}

    def __init__(self, seq):
        super().__init__(seq, SequenceType.RNA)

    def is_valid(self, seq):
        return all(i in RNASeq.nucleotides_dic for i in seq)

    # TODO: extract this method
    def gc_content(self):
        symbols = ['G', 'C']
        count = 0

        for i in self.seq:
            if i in symbols:
                count += 1
        return count/len(self.seq)

    # TODO: extract this method
    def reverse_complement(self):
        res = ''
        for base in reversed(self.seq):
            res += RNASeq.complement[base]
        
        return res

