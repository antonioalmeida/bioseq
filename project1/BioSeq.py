from SequenceType import SequenceType

class BioSeq:

    def __init__(self, seq, seq_type):
        self.seq = seq
        self.seq_type = seq_type

    def __str__(self):
        return self.seq + ' - ' + str(self.seq_type)
