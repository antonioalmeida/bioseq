import abc
import utils

TERMINATION = 0
DIAGONAL = 1
VERTICAL = 2
HORIZONTAL = 3

class Alignment():

    sm = []

    def __init__(self, sm=False):
        if sm:
            self.sm = sm
        else:
            self.sm = Alignment.read_substitution_matrix_file('bioseq/res/blosum62.mat')

    def run_global_align_multiple_solutions(self, seq1, seq2, g, debug=False):
        
        if debug:
            print('> Global Alignment between:')
            print('    > %s' % seq1)
            print('    > %s' % seq2)

        debug_2 = not (len(seq1) > 15 or len(seq2) > 15)
        (s,t) = self.global_align_multiple_solutions(seq1, seq2, g, debug_2)
        
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
    
        if debug:
            print('> Local Alignment between:')
            print('    > %s' % seq1)
            print('    > %s' % seq2)
            print()

        debug_2 = not (len(seq1) > 15 or len(seq2) > 15)
        (s,t,max_score) = self.local_align_multiple_solutions(seq1, seq2, g, debug_2)
        
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
        m = len(seq1)
        n = len(seq2)
        
        score = [[0]]
        trace = [[0]]
        
        # initialize gaps in rows
        score[0] = [g * j for j in range(0, n+1)]
        trace[0] = [[3] for _ in range(0, n+1)]
        
        # initialize gaps in cols
        for i in range(1, m+1):
            score.append([g*i])
            trace.append([[2]])

        # apply the recurrence to fill the matrices
        for i in range(1, m+1):
            for j in range(1, n+1):
                s1 = score[i-1][j-1] + self.__score_col_alignment(seq1[i-1], seq2[j-1], g)
                s2 = score[i-1][j] + g
                s3 = score[i][j-1] + g
                score[i].append(max(s1,s2,s3))
                trace[i].append(self.__max_indices(s1,s2,s3))
        
        if debug:
            print()
            print('> Score')
            utils.pretty_print_matrix(score)

            print()
            print('> Trace')
            utils.pretty_print_matrix(trace)
        
        return (score, trace)

    def recover_global_align_multiple_solutions(self, seq1, seq2, trace, g, debug=False):
        final = []
        alignments = [["","", len(seq1), len(seq2)]]

        while alignments:
            [al_1, al_2, i, j] = alignments.pop()
            if i == 0 & j == 0:
                final.append([al_1, al_2])
            else:
                for move in trace[i][j]:
                    next_al_1 = ""; next_al_2 = ""
                    next_i = i; next_j = j
                    if move is DIAGONAL:
                        next_al_1 = seq1[i-1] + al_1
                        next_al_2 = seq2[j-1] + al_2
                        next_i = i-1; next_j = j-1
                    elif move is HORIZONTAL:
                        next_al_1 = '_' + al_1
                        next_al_2 = seq2[j-1] + al_2
                        next_j = j-1
                    elif move is VERTICAL:
                        next_al_1 = seq1[i-1] + al_1
                        next_al_2 = '_' + al_2
                        next_i = i-1
                    new_alignment = [next_al_1, next_al_2, next_i, next_j]
                    alignments.append(new_alignment)

        return final

    def local_align_multiple_solutions(self, seq1, seq2, g, debug=False):
        """Local alignment"""
        score = [[0]]
        trace = [[[0]]]
        maxscore = 0
        # first row filled with zero
        for j in range(1, len(seq2)+1):
            score[0].append(0)
            trace[0].append([0])
        # first column filled with zero
        for i in range(1, len(seq1)+1):
            score.append([0])
            trace.append([[0]])
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = score[i][j] + self.__score_col_alignment(seq1[i], seq2[j], g)
                s2 = score[i][j+1] + g
                s3 = score[i+1][j] + g
                b = max(s1, s2, s3)
                if b <= 0:
                    score[i+1].append(0)
                    trace[i+1].append([0])
                else:
                    score[i+1].append(b)
                    trace[i+1].append(self.__max_indices(s1, s2, s3))
                    if b > maxscore:
                        maxscore = b

        if debug:
            print()
            print('> Score')
            utils.pretty_print_matrix(score)

            print()
            print('> Trace')
            utils.pretty_print_matrix(trace)
        return (score, trace, maxscore) 

    def recover_local_align_multiple_solutions(self, seq1, seq2, S, T, debug=False):
        final = []
        max_indices = self.__max_mat(S)
        alignments = [["","", i, j] for (i, j) in max_indices]

        if debug:
            print('> Number of maximum scores - %d' % len(max_indices))
            for (i,j) in max_indices:
                print('    > %d, %d' % (i,j))

        while alignments:
            [al_1, al_2, i, j] = alignments.pop()
            if T[i][j] == 0:
                final.append([al_1, al_2])
            else: 
                for move in T[i][j]:
                    next_al_1 = ""; next_al_2 = ""
                    next_i = i; next_j = j
                    if move is TERMINATION:
                        final.append([al_1, al_2])
                    else:
                        next_al_1 = ""; next_al_2 = ""
                        next_i = i; next_j = j
                        if move is DIAGONAL:
                            next_al_1 = seq1[i-1] + al_1
                            next_al_2 = seq2[j-1] + al_2
                            next_i = i-1; next_j = j-1
                        elif move == HORIZONTAL:
                            next_al_1 = '_' + al_1
                            next_al_2 = seq2[j-1] + al_2
                            next_j = j-1
                        elif move is VERTICAL:
                            next_al_1 = seq1[i-1] + al_1
                            next_al_2 = '_' + al_2
                            next_i = i-1
                        new_alignment = [next_al_1, next_al_2, next_i, next_j]
                        alignments.append(new_alignment)
        return final 

    def compare_pairwise_global_align(self, seqs, g):
        scores = [ [self.__get_global_alignment_max_score(seq1, seq2, g) for seq2 in seqs] for seq1 in seqs ]
        utils.pretty_print_with_header(scores, seqs)
        return scores

    def compare_pairwise_local_align(self, seqs, g):
        scores = [ [self.__get_local_alignment_max_score(seq1, seq2, g) for seq2 in seqs] for seq1 in seqs ]
        utils.pretty_print_with_header(scores, seqs)
        return scores

    def __max3t(self, v1,v2,v3):
        if v1 > v2:
            return 1 if v1 > v3 else 3
        else:
            return 2 if v2 > v3 else 3

    def __max_indices(self, *argv):
        m = max(argv)
        res = []
        for i,v in enumerate(argv):
            if(v == m):
                res.append(i+1)
        return res

    def __max_mat(self, mat):
        max_value = max([max(line) for line in mat])
        return [(i,j) for i,line in enumerate(mat) for j,v in enumerate(line) if v == max_value]

    # Provides the score of a column alignment, i.e., between characters c1 and c2
    # Assume a constant gap penalty g and a substitution matrix sm
    def __score_col_alignment(self, c1, c2, g):
        return g if c1 == '-' or c2 == '-' else self.sm[c1+c2]
    
    def __get_global_alignment_max_score(self, seq1, seq2, g):
        (s,_) = self.global_align_multiple_solutions(seq1, seq2, g)
        return s[len(seq1)][len(seq2)]

    def __get_local_alignment_max_score(self, seq1, seq2, g):
        (_,_,max_score) = self.local_align_multiple_solutions(seq1, seq2, g)
        return max_score

    @staticmethod
    def read_substitution_matrix_file(filename):
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
        dic = { x+y:match if x==y else mismatch for x in alphabet for y in alphabet}
        return dic
