from BioSeq import BioSeq
from RNASeq import RNASeq
from SequenceType import SequenceType

class DNASeq(BioSeq):

    nucleotides_dic = set(['A', 'T', 'C', 'G'])

    complement = {"A": "T",
                  "T": "A",
                  "G": "C",
                  "C": "G"}
    
    def __init__(self, seq):
        super().__init__(seq, SequenceType.DNA)

    def is_valid(self, seq):
        return all(i in DNASeq.nucleotides_dic for i in seq)

    def gc_content(self):
        symbols = ['G', 'C']
        count = 0

        for i in self.seq:
            if i in symbols:
                count += 1
        return count/len(self.seq)

    def reverse_complement(self):
        res = ''
        for base in reversed(self.seq):
            res += DNASeq.complement[base]
        
        return res
    
    def transcription(self):
        rna_seq = self.seq.replace('T','U')
        return RNASeq(rna_seq)

