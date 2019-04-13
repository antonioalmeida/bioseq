import abc
import utils

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

    def needleman_wunsch(self, seq1, seq2, g):
        m = len(seq1)
        n = len(seq2)
        
        score = [[0]]
        trace = [[0]]
        
        # initialize gaps in rows
        score[0] = [g * j for j in range(0, n+1)]
        trace[0] = [3 for _ in range(0, n+1)]
        
        
        # initialize gaps in cols
        for i in range(1, m+1):
            score.append([g*i])
            trace.append([2])

        # apply the recurrence to fill the matrices
        for i in range(1, m+1):
            for j in range(1, n+1):
                s1 = score[i-1][j-1] + self.__score_col_alignment(seq1[i-1], seq2[j-1], g)
                s2 = score[i-1][j] + g
                s3 = score[i][j-1] + g
                score[i].append(max(s1,s2,s3))
                trace[i].append(self.__max3t(s1,s2,s3))
        
        return (score, trace)

    def recover_global_align(self, trace, seq1, seq2):
        res = ["",""]
        (i,j) = (len(seq1),len(seq2))
        
        while i>0 or j>0:
            if trace[i][j] is DIAGONAL:
                res[0] += seq1[i-1]
                res[1] += seq2[j-1]
                i -= 1; j -= 1
            elif trace[i][j] is HORIZONTAL:
                res[0] += '_'
                res[1] += seq2[j-1]
                j -= 1
            elif trace[i][j] is VERTICAL:
                res[0] += seq1[i-1]
                res[1] += '_'
                i -= 1
        
        return [s[::-1] for s in res]

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

    # Provides the score of a column alignment, i.e., between characters c1 and c2
    # Assume a constant gap penalty g and a substitution matrix sm
    def __score_col_alignment(self, c1, c2, g):
        return g if c1 == '-' or c2 == '-' else self.sm[c1+c2]

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
