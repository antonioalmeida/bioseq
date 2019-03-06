from bioseq.bioseq import BioSeq, SequenceType
from bioseq.dnaseq import DNASeq
from bioseq.rnaseq import RNASeq
from bioseq.proteinseq import AminoacidSeq, ProteinSeq

__types_dic = {
    'dna': DNASeq,
    'rna': RNASeq,
    'aa': AminoacidSeq,
    'aminoacid': AminoacidSeq,
    'protein': ProteinSeq,
}

def read_fasta_file(filename, sequence_type='dna'):
    dic = __create_dictionary(filename)

    # TODO: verify sequence_type input
    Type = __types_dic[sequence_type]
    return {k:Type(v) for k,v in dic.items()} 

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