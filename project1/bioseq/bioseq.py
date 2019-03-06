import abc
from enum import Enum 

class SequenceType(Enum):
    DNA = 'DNA'
    RNA = 'RNA'
    AMINOACID = 'AMINOACID'
    PROTEIN = 'PROTEIN'
    
    def __str__(self):
        return format(self.value)
    
class BioSeq(abc.ABC):
    """Abstract class representing a biological sequence."""

    def __init__(self, seq, seq_type):
        assert self.is_valid(seq), 'ERROR: invalid sequence'
        self.seq = seq.upper()
        self.seq_type = seq_type

    @abc.abstractmethod
    def is_valid(self, seq):
        """Checks if the sequence is valid, depending on its type."""

    def symbol_frequency(self):
        """Calculates the symbol frequency of the sequence."""
        dic = {}
        for i in self.seq:
            if(i not in dic):
                dic[i] = 1
            else:
                dic[i] += 1
        return dic
    
    def pretty_print(self):
        """Prints to the console summarized information about the sequence."""

        sf = self.symbol_frequency()
        most_common = max(sf, key=(lambda x : sf[x]))

        print("")
        print("Sequence")
        print("    Type: " + str(self.seq_type))
        print("    Length: " + str(len(self.seq)) + " characters")
        print("    Symbols: contains " + str(len(sf)) + " different symbols")
        print("           : the most common symbol is " + most_common)
        print("    First 100 characters: " + self.seq[0:100])
        print("")
    
    def save_to_file(self, filename):
        import pickle
        with open(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
            print('Saved sequence object to ' + filename)

    def __str__(self):
        return self.seq + ' - ' + str(self.seq_type)

    def __eq__(self, other):
        return self.seq == other.seq and self.seq_type == other.seq_type

    def __len__(self):
        return len(self.seq)
