from enum import Enum 

class SequenceType(Enum):
    DNA = 'DNA'
    RNA = 'RNA'
    PROTEIN = 'PROTEIN'
    
    def __str__(self):
        return format(self.value)
    