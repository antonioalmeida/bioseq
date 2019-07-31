from bioseq.bioseq import BioSeq, SequenceType
from bioseq.dnaseq import DNASeq
from bioseq.rnaseq import RNASeq
from bioseq.proteinseq import AminoacidSeq, ProteinSeq
import re

__types_dic = {
    'dna': DNASeq,
    'rna': RNASeq,
    'aa': AminoacidSeq,
    'aminoacid': AminoacidSeq,
    'protein': ProteinSeq,
}

def read_fasta_file(filename, sequence_type='dna', use_dic=True):
    dic = __create_dictionary(filename)
    Type = __types_dic[sequence_type]

    if use_dic:
        return {k:Type(v) for k,v in dic.items()} 
    else:
        return [Type(v, k) for k,v in dic.items()]

def __create_dictionary(filename):
    dic = {}
    fd = open(filename)
    key = ''
    
    for line in fd:
        if len(line.strip()):
            line = line.strip()

            if(line[0] == '>'):
                key = line[1:]
                dic[key] = ''
            else:
                dic[key] += line
    
    fd.close()
    return dic