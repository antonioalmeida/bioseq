import abc

class Alignment():

    sm = []

    def __init__(self, sm_file='bioseq/res/blosum62.mat'):
        self.sm = self.__read_substitution_matrix_file(sm_file)

    def needleman_wunsch(self, seq1, seq2, sm, g):
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
                s1 = score[i-1][j-1] + self.__score_col_alignment(seq1[i-1], seq2[j-1], sm, g)
                s2 = score[i-1][j] + g
                s3 = score[i][j-1] + g
                score[i].append(max(s1,s2,s3))
                trace[i].append(self.__max3t(s1,s2,s3))
        
        return (score, trace)

    def __max3t(self, v1,v2,v3):
        if v1 > v2:
            return 1 if v1 > v3 else 3
        else:
            return 2 if v2 > v3 else 3

    # Provides the score of a column alignment, i.e., between characters c1 and c2
    # Assume a constant gap penalty g and a substitution matrix sm
    def __score_col_alignment(self, c1, c2, sm, g):
        return g if c1 == '-' or c2 == '-' else sm[c1+c2]

    
    def __read_substitution_matrix_file(self, filename):
        
        with open(filename) as f:
            
            chars = f.readline().strip().replace('\t','')
        
            keys = [x for x in chars]
            matrix = []
            
            for line in f: 
                matrix.append(line.strip().split())
            
            dic = { x+y:int(matrix[yi][xi]) for xi,x in enumerate(keys) for yi,y in enumerate(keys) }
            return dic