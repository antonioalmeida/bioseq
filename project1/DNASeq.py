from BioSeq import BioSeq
from SequenceType import SequenceType

class DNASeq(BioSeq):
    
    def __init__(self, seq):
        super().__init__(seq, SequenceType.DNA)
