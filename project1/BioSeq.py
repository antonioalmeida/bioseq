import abc
from SequenceType import SequenceType

class BioSeq(abc.ABC):

    # Constructor 
    def __init__(self, seq, seq_type):
        self.seq = seq.upper()
        self.seq_type = seq_type

    @abc.abstractmethod
    def is_valid(self):
        """Checks if the sequence is valid, depending on its type."""

    def symbol_frequency(self):
        """Frequency of Symbols"""
        dic = {}
        for i in self.seq:
            if(i not in dic):
                dic[i] = 1
            else:
                dic[i] += 1
        return dic

    @abc.abstractmethod
    def gc_content(self):
        """Calculate the frequency of G and C nucleotides (percentage of 'G' and 'C') of the sequence."""

    def __str__(self):
        return self.seq + ' - ' + str(self.seq_type)

