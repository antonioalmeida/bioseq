from BioSeq import BioSeq
from SequenceType import SequenceType

class ProteinSeq(BioSeq):

    def __init__(self, seq):
        super().__init__(seq, SequenceType.PROTEIN)

    def is_valid(self, seq):
        return all(i.isalpha() or i in ['_', '*'] for i in seq)
