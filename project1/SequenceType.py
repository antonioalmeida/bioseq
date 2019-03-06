from enum import Enum 

class SequenceType(Enum):
    DNA = 'DNA'
    RNA = 'RNA'
    AMINOACID = 'AMINOACID'
    PROTEIN = 'PROTEIN'
    
    def __str__(self):
        return format(self.value)
    