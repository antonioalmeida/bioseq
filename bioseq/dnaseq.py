from bioseq.bioseq import SequenceType
from bioseq.dnracommon import DRNACommon
from bioseq.rnaseq import RNASeq

class DNASeq(DRNACommon):

    aminoacids_dic = DRNACommon.read_dic_aminoacids('bioseq/res/CODON_AMINOACID_MAP_DNA.txt')

    nucleotides_dic = set(['A', 'T', 'C', 'G'])

    complement = {"A": "T",
                  "T": "A",
                  "G": "C",
                  "C": "G"}

    def __init__(self, seq):
        super().__init__(seq, SequenceType.DNA)

    def transcription(self):
        """Computes the correspondent RNA sequence corresponding to the transcription of the DNA sequence.
        
        Returns:
            An instance of RNASeq.
        """
        rna_seq = self.seq.replace('T','U')
        return RNASeq(rna_seq)
