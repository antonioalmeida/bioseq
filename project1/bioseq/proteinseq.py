from bioseq.bioseq import BioSeq, SequenceType

class AminoacidSeq(BioSeq):

    def __init__(self, seq):
        super().__init__(seq, SequenceType.AMINOACID)

    def is_valid(self, seq):
        return all(i.isalpha() or i in ['_', '*'] for i in seq)

    def all_proteins_rf(self, min_size=0):
        """Computes all possible proteins in an aminoacid sequence. 

        Args:
            min_size: Minimum length of each protein sequence.
        
        Returns:
            A list of instances of ProteinSeq, representing the possible proteins. Sorted by decreasing sequence length.
        """
        proteins = []
        index_stack = []
    
        for i,x in enumerate(self.seq):
            if x == 'M': # start_sequence
                index_stack.append(i)
            if x == '_':
                for index in index_stack:
                    proteins.append(ProteinSeq(self.seq[index:i+1]))
                index_stack = []
    
        proteins = [p for p in proteins if len(p) > min_size]
        proteins = sorted(proteins, key=len, reverse=True)
        return proteins

class ProteinSeq(BioSeq):

    def __init__(self, seq):
        super().__init__(seq, SequenceType.PROTEIN)

    def is_valid(self, seq):
        return all(i.isalpha() or i in ['_', '*'] for i in seq) and seq[0] == 'M' and seq[len(seq) - 1] == '_'
