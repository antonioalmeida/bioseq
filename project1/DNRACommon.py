import abc
from BioSeq import BioSeq
from ProteinSeq import ProteinSeq
from SequenceType import SequenceType

class DRNACommon(BioSeq):

    aminoacids_dic = {}
    nucleotides_dic = {}
    complement = {}

    def __init__(self, seq, seq_type):
        super().__init__(seq, seq_type)

    def is_valid(self, seq):
        return all(i in self.nucleotides_dic for i in seq)

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
            res += self.complement[base]
        
        return res

    def translation(self, init_index=0):
        res = ''

        for i in range(init_index, len(self.seq) - 2, 3):
            res += self.__translate_codon(self.seq[i:i+3])
        return ProteinSeq(res)
    
    def __translate_codon(self, triplet):
        assert triplet in self.aminoacids_dic, 'Triplet not found in aminoacids dictionary'
        return self.aminoacids_dic[triplet]

    @staticmethod
    def read_dic_aminoacids(filename):
        dic = {}
        fd = open(filename)
        
        for line in fd:
            line = line.strip()
            triplet = line[1:4]
            aminoacid = line[7]
            dic[triplet] = aminoacid
        
        return dic



