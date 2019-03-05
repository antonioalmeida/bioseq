from BioSeq import BioSeq
from DNRACommon import DRNACommon
from SequenceType import SequenceType

class RNASeq(DRNACommon):

    aminoacids_dic = DRNACommon.read_dic_aminoacids('CODON_AMINOACID_MAP_RNA.txt')

    nucleotides_dic = set(['A', 'U', 'C', 'G'])

    complement = {"A": "U",
                  "U": "A",
                  "G": "C",
                  "C": "G"}

    def __init__(self, seq):
        super().__init__(seq, SequenceType.RNA)

    def is_valid(self, seq):
        return all(i in RNASeq.nucleotides_dic for i in seq)

    # TODO: this is wrong
    def transcription(self):
        rna_seq = self.seq.replace('T','U')
        return RNASeq(rna_seq)

