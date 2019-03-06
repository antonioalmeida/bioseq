import abc
from bioseq.bioseq import BioSeq, SequenceType
from bioseq.proteinseq import AminoacidSeq

class DRNACommon(BioSeq):
    "Abstract class to encapsulate the common methods of DNASeq and RNASeq."

    aminoacids_dic = {}
    nucleotides_dic = {}
    complement = {}

    def __init__(self, seq, seq_type):
        super().__init__(seq, seq_type)

    def is_valid(self, seq):
        return all(i in self.nucleotides_dic for i in seq)

    def gc_content(self):
        """Calculate the frequency of G and C nucleotides (percentage of 'G' and 'C') of the sequence."""
        symbols = ['G', 'C']
        count = 0

        for i in self.seq:
            if i in symbols:
                count += 1
        return count/len(self.seq)

    def reverse_complement(self):
        """Calculates the reverse complement of a DNA molecule."""
        res = ''
        for base in reversed(self.seq):
            res += self.complement[base]

        return res

    def translation(self, init_index=0):
        """Translates a DNA or RNA sequence into an aminoacid sequence.

        Args:
            init_index: Index in the sequence where the translation starts. Defaults to 0.
        Returns:
            The translated sequence, as an instance of AminoacidSeq."""
        res = self.__translation_helper(self.seq, init_index)
        return AminoacidSeq(res)

    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid, in a DNA  or RNA sequence.
        
        Args:
            aminoacid: A string character representing the target aminoacid.
        Returns:
            A dictionary with the frequency of each codon.
        """
        res_dic = {k:0 for k,v in self.aminoacids_dic.items() if v is aminoacid}

        for i in range(0, len(self.seq) - 2, 3):
            codon = self.seq[i:i+3]
            if self.__translate_codon(codon) == aminoacid:
                res_dic[codon] += 1

        return res_dic


    def reading_frames(self):
        """Computes the six reading frames of a DNA or RNA sequence and its reverse complement.
        
        Returns:
            A list of instances of AminoacidSeq representing the sequence's reading frames."""
        res = []
        for i in range(0,3):
            res.append(self.translation(i))
        
        rev_comp = self.reverse_complement()
    
        for i in range(0,3):
            rf = self.__translation_helper(rev_comp, i)
            res.append(AminoacidSeq(rf))
    
        return res

    def all_open_reading_frames(self, min_size=0):
        """Computes all putative proteins of an DNA or RNA sequence.
        
        Args:
            min_size: Minimum length of each protein sequence.
        
        Returns:
            A list of instances of ProteinSeq, representing the possible proteins. Sorted by decreasing sequence length.
        """
        rfs = self.reading_frames()
        possible_proteins = []
    
        for rf in rfs:
            possible_proteins += rf.all_proteins_rf(min_size)
        
        sorted_proteins = sorted(possible_proteins, key=len, reverse=True)
        return sorted_proteins

    def __translate_codon(self, triplet):
        """Translates a codon into an aminoacid using the standard genetic code dictionary.

        Args:
            triplet: A character representing the target triplet.

        Returns:
            The correspondent aminoacid.
        """
        assert triplet in self.aminoacids_dic, 'Triplet not found in aminoacids dictionary'
        return self.aminoacids_dic[triplet]

    def __translation_helper(self, seq, init_index=0):
        res = ''
        for i in range(init_index, len(seq) - 2, 3):
            res += self.__translate_codon(seq[i:i+3])
        return res

    @staticmethod
    def read_dic_aminoacids(filename):
        dic = {}
        fd = open(filename)

        for line in fd:
            line = line.strip()
            triplet = line[1:4]
            aminoacid = line[7]
            dic[triplet] = aminoacid

        fd.close()
        return dic
