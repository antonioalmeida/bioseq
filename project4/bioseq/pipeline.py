from bioseq.fasta import read_fasta_file
import bioseq.blast
import bioseq.msa
import bioseq.upgma

class Pipeline():

    db_filename = ''
    seq = ''
    w = 0

    def __init__(self, query_filename, type='dna', db_filename='bioseq/res/seqdump.txt', w=3, seq_type='protein'):
        self.seq = read_fasta_file(query_filename, seq_type, use_dic=False)[0]
        self.db_filename = db_filename
        self.w = w

    def run_blast(self):
        b = bioseq.blast(self.db_filename, self.w)
        r = b.best_n_alignments(self.seq)
        return r

    def run_msa(self, seqs, sm=False, g=-3):
        msa = bioseq.msa(seqs, sm, g, species=True)
        (consensus, alignment) = msa.align()

        score = msa.score_sp(alignment)
        alignment.pretty_print(score)
        return alignment

    def run_upgma(self, seqs):
        t = bioseq.upgma(seqs).execute()
        #t.print_tree_2()
        #//t.print_tree()
        print(t)





