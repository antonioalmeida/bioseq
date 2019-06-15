from bioseq.fasta import read_fasta_file
import bioseq.blast

class Pipeline():

    db_filename = ''
    desc = ''
    seq = ''
    w = 0

    def __init__(self, query_filename, type='dna', db_filename='bioseq/res/seqdump.txt', w=3, seq_type='protein'):
        self.desc, self.seq = list(read_fasta_file(query_filename, seq_type).items())[0]
        self.db_filename = db_filename
        self.w = w

    def run_blast(self):
        b = bioseq.blast(self.db_filename, self.w)
        r = b.best_n_alignments(self.seq)
        return r


        

