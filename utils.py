def pretty_print_with_header(m, seqs):
    matrix = []
    header = [' ']
    header += (seqs)
    matrix.append(header)

    for i,line in enumerate(m):
        n_line = [seqs[i][:3]]
        n_line += line
        matrix.append(n_line)

    pretty_print_matrix(matrix)

def pretty_print_matrix(m):
    s = [[str(e) for e in row] for row in m]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))

def print_tc(score, trace):
    print()
    print('> Score')
    pretty_print_matrix(score)

    print()
    print('> Trace')
    pretty_print_matrix(trace)
