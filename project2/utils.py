
def pretty_print(seq1, seq2, m):
    
    for line in m:
        pretty_row(line)

def pretty_header(seq):
    print('    ', end="")
    for i in seq:
        print(i + ' | ', end="")
    print()

def pretty_row(row):
    for i in row:
        print(str(i) + ' | ', end="")
    print()

def pretty_print_matrix(m):
    s = [[str(e) for e in row] for row in m]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))