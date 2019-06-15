from bioseq.fasta import read_fasta_file

class Pipeline():

    desc = ""
    seq = ""

    def __init__(self, filename, seq_type='dna'):
        self.desc, self.seq = list(read_fasta_file(filename, seq_type).items())[0]


        

