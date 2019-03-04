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

    def translation(self, init_index=0):
        res = ''

        for i in range(init_index, len(self.seq) - 2, 3):
            res += self.__translate_codon(self.seq[i:i+3])
        return res

    def __translate_codon(self, triplet):
        assert triplet in DNASeq.__aminoacids_dic, 'Triplet not found in aminoacids dictionary'
        return DNASeq.__aminoacids_dic[triplet]

    def read_dic_aminoacids(filename):
        dic = {}
        fd = open(filename)
        
        for line in fd:
            line = line.strip()
            triplet = line[1:4]
            aminoacid = line[7]
            dic[triplet] = aminoacid
        
        return dic

    __aminoacids_dic = read_dic_aminoacids('genetic_code.txt')



