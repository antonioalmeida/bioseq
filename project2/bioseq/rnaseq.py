from bioseq.bioseq import BioSeq, SequenceType
from bioseq.dnracommon import DRNACommon
import bioseq.dnaseq

class RNASeq(DRNACommon):

    aminoacids_dic = DRNACommon.read_dic_aminoacids('bioseq/res/CODON_AMINOACID_MAP_RNA.txt')

    nucleotides_dic = set(['A', 'U', 'C', 'G'])

    complement = {"A": "U",
                  "U": "A",
                  "G": "C",
                  "C": "G"}

    def __init__(self, seq):
        super().__init__(seq, SequenceType.RNA)

    def is_valid(self, seq):
        return all(i in RNASeq.nucleotides_dic for i in seq)

    def reverse_transcription(self):
        """Computes the correspondent DNA sequence corresponding to the reverse transcription of the RNA sequence.
        
        Returns:
            An instance of DNASeq.
        """
        dna_seq = self.seq.replace('U','T')
        return bioseq.dnaseq.DNASeq(dna_seq)
