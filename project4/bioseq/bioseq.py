import abc
from enum import Enum 
import pickle
import re

class SequenceType(Enum):
    DNA = 'DNA'
    RNA = 'RNA'
    AMINOACID = 'AMINOACID'
    PROTEIN = 'PROTEIN'
    
    def __str__(self):
        return format(self.value)
    
class BioSeq(abc.ABC):
    """Abstract class representing a biological sequence."""

    def __init__(self, seq, seq_type, desc=''):
        assert self.is_valid(seq), 'ERROR: invalid sequence'
        self.seq = seq.upper()
        self.seq_type = seq_type
        self.desc = desc

        m = re.search(r"\[(.*?)\]", desc)
        if m:
            self.species = m.group(0)
        else:
            self.species = ''

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
        print("    Description: " + str(self.desc))
        print("    Species: " + str(self.species))
        print("    Length: " + str(len(self.seq)) + " characters")
        print("    Symbols: contains " + str(len(sf)) + " different symbols")
        print("           : the most common symbol is " + most_common)
        print("    First 100 characters: " + self.seq[0:100])
        print("")
    
    def save_to_file(self, filename):
        with open(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
            print('Saved sequence object to ' + filename)

    @staticmethod
    def load_from_file(filename):
        with open(filename, 'rb') as input:
            print('Loading sequence object from ' + filename + '...')
            return pickle.load(input)

    def __str__(self):
        return self.species if self.species != '' else self.seq + ' - ' + str(self.seq_type)

    def __eq__(self, other):
        return self.seq == other.seq and self.seq_type == other.seq_type

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        return self.seq.__getitem__(key)
