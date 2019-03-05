from BioSeq import BioSeq
from RNASeq import RNASeq
from DNRACommon import DRNACommon
from SequenceType import SequenceType

class DNASeq(DRNACommon):

    aminoacids_dic = DRNACommon.read_dic_aminoacids('CODON_AMINOACID_MAP_DNA.txt')

    nucleotides_dic = set(['A', 'T', 'C', 'G'])

    complement = {"A": "T",
                  "T": "A",
                  "G": "C",
                  "C": "G"}

    def __init__(self, seq):
        super().__init__(seq, SequenceType.DNA)

    def transcription(self):
        rna_seq = self.seq.replace('T','U')
        return RNASeq(rna_seq)