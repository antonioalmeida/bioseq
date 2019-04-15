import abc
import utils

GAP = '_'
TERMINATION = 0
DIAGONAL = 1
VERTICAL = 2
HORIZONTAL = 3

class Alignment():
    sm = [] # each Alignment issue has an implicit substitution matrix

    def __init__(self, sm=False):
        """Constructor for the Alignment class"""
        if sm:
            self.sm = sm
        else: # uses blosum62 as sm by default
            self.sm = Alignment.read_substitution_matrix_file('bioseq/res/blosum62.mat')

    def run_global_align_multiple_solutions(self, seq1, seq2, g, debug=False):
        """Runs and recovers the global alignment solutions"""
        if debug:
            print('> Global Alignment between:')
            print('    > %s' % seq1)
            print('    > %s' % seq2)

        debug_2 = debug and not (len(seq1) > 15 or len(seq2) > 15)
        (_,t) = self.global_align_multiple_solutions(seq1, seq2, g, debug_2)
        alignments = self.recover_global_align_multiple_solutions(seq1, seq2, t, g, debug)

        if debug:
            print('> Got %d alignment(s)' % len(alignments))
            print('> Some results:')
            for [al1, al2] in alignments[:4]:
                print('    > %s' % al1)
                print('    > %s' % al2)
                print()
        
        return alignments

    def run_local_align_multiple_solutions(self, seq1, seq2, g, debug=False):
        """Runs and recovers the local alignment solutions"""
        if debug:
            print('> Local Alignment between:')
            print('    > %s' % seq1)
            print('    > %s' % seq2)
            print()

        debug_2 = debug and not (len(seq1) > 15 or len(seq2) > 15)
        (s,t,_) = self.local_align_multiple_solutions(seq1, seq2, g, debug_2)
        alignments = self.recover_local_align_multiple_solutions(seq1, seq2, s, t, debug)

        if debug:
            print('> Got %d alignments' % len(alignments))
            print('> Some results:')
            for [al1, al2] in alignments[:4]:
                print('    > %s' % al1)
                print('    > %s' % al2)
                print()
        
        return alignments
        
    def global_align_multiple_solutions(self, seq1, seq2, g, debug=False):
        """Runs the Needleman–Wunsch algorithm between seq1 and seq2, with a gap of g"""
        m = len(seq1); n = len(seq2)
        S = [[0]]; T = [[0]]
        
        # initialize gaps in rows
        S[0] = [g * j for j in range(0, n+1)]
        T[0] = [[3] for _ in range(0, n+1)]
        
        # initialize gaps in cols
        for i in range(1, m+1):
            S.append([g*i])
            T.append([[2]])

        # apply the recurrence to fill the matrices
        for i in range(1, m+1):
            for j in range(1, n+1):
                s1 = S[i-1][j-1] + self._score_col_alignment(seq1[i-1], seq2[j-1], g)
                s2 = S[i-1][j] + g
                s3 = S[i][j-1] + g
                S[i].append(max(s1,s2,s3))
                T[i].append(self._max_indices(s1,s2,s3))
        
        if debug: utils.print_tc(S, T)
        return (S, T)

    def recover_global_align_multiple_solutions(self, seq1, seq2, trace, g, debug=False):
        """Runs the traceback algorithm to recover the solution for global alignment between seq1 and seq2"""
        final = []
        alignments = [["","", len(seq1), len(seq2)]]

        while alignments:
            [al_1, al_2, i, j] = alignments.pop()
            if i == 0 & j == 0:
                final.append([al_1, al_2])
            else:
                for move in trace[i][j]:
                    alignments.append(self._recover_step(move, seq1, seq2, al_1, al_2, i, j))
        return final

    def local_align_multiple_solutions(self, seq1, seq2, g, debug=False):
        """Runs the Smith-Waterman algorithm between seq1 and seq2, with a gap of g"""
        m = len(seq1); n = len(seq2)
        S = [[0]]; T = [[[0]]]
        max_s = 0

        for j in range(1, n+1):
            S[0].append(0)
            T[0].append([0])

        for i in range(1, m+1):
            S.append([0])
            T.append([[0]])

        for i in range(0, m):
            for j in range(n):
                s1 = S[i][j] + self._score_col_alignment(seq1[i], seq2[j], g)
                s2 = S[i][j+1] + g
                s3 = S[i+1][j] + g
                b = max(s1, s2, s3)
                if b <= TERMINATION:
                    S[i+1].append(0)
                    T[i+1].append([0])
                else:
                    S[i+1].append(b)
                    T[i+1].append(self._max_indices(s1, s2, s3))
                    if b > max_s:
                        max_s = b

        if debug: utils.print_tc(S, T)
        return (S, T, max_s) 

    def recover_local_align_multiple_solutions(self, seq1, seq2, S, T, debug=False):
        """Runs the traceback algorithm to recover the solution for local alignment between seq1 and seq2"""
        final = []
        max_indices = self._max_mat(S)
        alignments = [["","", i, j] for (i, j) in max_indices]

        if debug:
            print('> Number of maximum scores - %d' % len(max_indices))
            for (i,j) in max_indices:
                print('    > %d, %d' % (i,j))

        while alignments:
            [al_1, al_2, i, j] = alignments.pop()
            for move in T[i][j]:
                if move == TERMINATION:
                    final.append([al_1, al_2])
                else:
                    alignments.append(self._recover_step(move, seq1, seq2, al_1, al_2, i, j))
        return final 

    def compare_pairwise_global_align(self, seqs, g):
        """Creates a matrix with the global alignment scores between all sequences on seqs, using gap g"""
        scores = [ [self._get_global_alignment_max_score(seq1, seq2, g) for seq2 in seqs] for seq1 in seqs ]
        utils.pretty_print_with_header(scores, seqs)
        return scores

    def compare_pairwise_local_align(self, seqs, g):
        """Creates a matrix with the local alignment scores between all sequences on seqs, using gap g"""
        scores = [ [self._get_local_alignment_max_score(seq1, seq2, g) for seq2 in seqs] for seq1 in seqs ]
        utils.pretty_print_with_header(scores, seqs)
        return scores

    def _recover_step(self, move, s1, s2, l, r, i, j):
        """Helper function to isolate the recover step on local and global alignments"""
        if move == DIAGONAL:
            l = s1[i-1] + l; r = s2[j-1] + r
            i -= 1; j-= 1
        elif move == HORIZONTAL:
            l = GAP + l; r = s2[j-1] + r
            j-= 1
        elif move == VERTICAL:
            l = s1[i-1] + l; r = GAP + r
            i-= 1
        return (l,r,i,j)

    def _max_indices(self, *argv):
        """Returns the indices of the max values on argv"""
        m = max(argv)
        res = []
        for i,v in enumerate(argv):
            if(v == m):
                res.append(i+1)
        return res

    def _max_mat(self, mat):
        """Returns indices on matrix mat that hold the greatest value"""
        max_value = max([max(line) for line in mat])
        return [(i,j) for i,line in enumerate(mat) for j,v in enumerate(line) if v == max_value]

    def _score_col_alignment(self, c1, c2, g):
        """Provides the score of a column alignment, i.e., between characters c1 and c2"""
        return g if c1 == '-' or c2 == '-' else self.sm[c1+c2]
    
    def _get_global_alignment_max_score(self, seq1, seq2, g):
        """Run the global alignment algorithm and returns the max score"""
        (s,_) = self.global_align_multiple_solutions(seq1, seq2, g)
        return s[len(seq1)][len(seq2)]

    def _get_local_alignment_max_score(self, seq1, seq2, g):
        """Runs the local alignment algorithm and returns the max score"""
        (_,_,max_score) = self.local_align_multiple_solutions(seq1, seq2, g)
        return max_score

    @staticmethod
    def read_substitution_matrix_file(filename):
        """Reads substitutions matrix from filename"""
        with open(filename) as f:
            chars = f.readline().strip().replace('\t','')
            keys = [x for x in chars]
            matrix = []
            for line in f: 
                matrix.append(line.strip().split())
            dic = { x+y:int(matrix[yi][xi]) for xi,x in enumerate(keys) for yi,y in enumerate(keys) }
            return dic

    @staticmethod 
    def create_substitution_matrix(alphabet, match=1, mismatch=0):
        """Creates substituion matrix based on an alphabet, a match value, and a mismatch value"""
        dic = { x+y:match if x==y else mismatch for x in alphabet for y in alphabet}
        return dic
