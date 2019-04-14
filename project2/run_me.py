import bioseq
import itertools

def press_enter():
    print()
    input('Press ENTER to continue...')
    for i in range(0,50): print()

print('Welcome to the bioseq v2 package demo.')

print('> C1')
print('We\'ll now be calculating global and local alignments from protein_sequences.fas file')

# Test C1
al = bioseq.alignment() # this uses bosum62 by default
fasta_dic = bioseq.read_fasta_file('examples/protein_sequences.fas', 'aminoacid') 
pairs = list(itertools.combinations(fasta_dic.keys(), 2))
alignments = []

print('Calculating alignments...')
for (key1, key2) in pairs:
    seq1 = fasta_dic[key1]; seq2 = fasta_dic[key2]
    n_global_alignments = len(al.run_global_align_multiple_solutions(seq1, seq2, -3))

    (s,t, _) = al.local_align_multiple_solutions(seq1, seq2, -3)
    n_local_alignments = len(al.recover_local_align_multiple_solutions(seq1, seq2,s,t))

    alignments.append([key1, key2, n_global_alignments, n_local_alignments])

sorted_global = sorted(alignments, key=lambda x: x[2], reverse=True)
sorted_local = sorted(alignments, key=lambda x: x[3], reverse=True)

print('> Top 5 global alignments on protein_sequences.fas')
for line in sorted_global[:5]:
    print('    > %s %s - %d' % (line[0], line[1], line[2]))

print('> Top 5 local alignments on protein_sequences.fas')
for line in sorted_local[:5]:
    print('    > %s %s - %d' % (line[0], line[1], line[3]))
press_enter()

seq1 = 'GATTACA'
seq2 = 'GCATGCT'
print('> C2')
print('We\'ll now be calculating local and global alignments for %s and %s, with match = 1, mismatch = -1, and gap = -1' % (seq1, seq2))
press_enter()

print('First, global alignment...')
press_enter()
sm = bioseq.alignment.create_substitution_matrix('ATCG', 1, -1)
al = bioseq.alignment(sm)
al.run_global_align_multiple_solutions(seq1, seq2, -1, debug=True)
press_enter()

print('Now, local alignment...')
press_enter()
al.run_local_align_multiple_solutions(seq1, seq2, -1, debug=True)
press_enter()

print('> D')
print('We\'ll now be showcasing compare_pairwise_global_align with a few sample sequences')
press_enter()
al = bioseq.alignment()
seqs = ['AGT', 'ACG','TAG', 'CCT']
al.compare_pairwise_global_align(seqs, -3)
press_enter()

print('We\'ll now be showcasing compare_pairwise_local_align with a few sample sequences')
press_enter()
al = bioseq.alignment()
seqs = ['AGT', 'ACG','TAG', 'CCT']
al.compare_pairwise_local_align(seqs, -3)
press_enter()