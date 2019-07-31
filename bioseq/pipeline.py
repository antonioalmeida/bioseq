from bioseq.fasta import read_fasta_file
import bioseq.blast
import bioseq.msa
import bioseq.upgma
import bioseq.graph

class Pipeline():

    db_filename = ''
    seq = ''
    w = 0

    blast = None
    tree = None
    msa = None 
    graph = None
    distances = None

    def __init__(self, query_filename, type='dna', db_filename='bioseq/res/seqdump.txt', w=3, seq_type='protein'):
        self.seq = read_fasta_file(query_filename, seq_type, use_dic=False)[0]
        self.db_filename = db_filename
        self.w = w

    def run_blast(self):
        b = bioseq.blast(self.db_filename, self.w)
        print('> Executing BLAST...')
        seqs, scores = b.best_n_alignments(self.seq)
        self.blast = seqs
        print('> BLAST execution finished.')
        
        for i in range(len(scores)):
            print('   > (%s, score=%d)' % (str(seqs[i]), scores[i]))
        return seqs

    def run_msa(self, seqs=None, sm=False, g=-3):
        if not seqs and not self.blast:
            print('> Execute <obj>.run_blast() first.'); return
        if not seqs:
            seqs = self.blast

        print('> Executing MSA...')
        msa = bioseq.msa(seqs + [self.seq], sm, g, species=True)
        (consensus, alignment) = msa.align()
        self.msa = alignment

        print('> MSA execution finished.')
        print('    > You can now use <obj>.visualize_msa() to view the resulting alignment.')
        return alignment

    def run_upgma(self, seqs=None):
        if not seqs and not self.msa:
            print('> Execute <obj>.run_msa() first.'); return
        if not seqs:
            seqs = self.msa

        up = bioseq.upgma(seqs)
        self.distances = up.original_dist
        self.tree = up.execute()
        print('> UPGMA execution finished.')
        print('     > You can now use <obj>.phylo() and <obj>.text() to visualize the tree.')

    def run_network(self, cut=15):
        if self.distances: 
            self.network = bioseq.Graph(m=self.distances, c=cut)
        else: 
            print('> Execute <obj>.run_upgma() first.')
            return None

        print('> Network creation finished.')
        print('     > You can now use <obj>.stats() and <obj>.render() to interact with the network.')
        return self.network

    def export_network(self, cut=15):
        if self.distances: 
            network = bioseq.Graph(m=self.distances, c=cut)
            path='out/cut-%d' % cut
            network.export(path)
            print('> Graph exported to %s' % path)
        else: 
            print('> Execute <obj>.run_upgma() first.')
            return None

    def visualize_msa(self):
        if self.msa: self.msa.print()
        else: print('> Execute <obj>.run_msa() first.')

    def phylo(self):
        if(self.tree): self.tree.phylo_tree()
        else: print('> Execute <obj>.run_upgma() first.')

    def text(self):
        if(self.tree): self.tree.print_tree()
        else: print('> Execute <obj>.run_upgma() first.')








