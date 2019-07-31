import bioseq as bs
import itertools

def press_enter():
    print()
    input('Press ENTER to continue...')
    for i in range(0,100): print()

print('Welcome to the bioseq v3 package demo.')

print('> BLAST')
print('     > We\'ll now be demonstrating BLAST execution on the given database.')
print('     > Using blosum62 and w=3')
press_enter()

p = bs.pipeline('bioseq/res/source.fasta')
p.run_blast()
press_enter()

print('> MSA')
print('     > We\'ll now be executing MSA between the resulting sequences and the query sequence. (Homo Sapiens)')
print('     > Using blosum62 and g=-3')
press_enter()

p.run_msa()

print()
print('     > This alignment is hard to visualize, so instead, here\'s a simpler MSA visualization.')
print('     > MSA between "ATAGC", "AACC", "ATGAC", simple match substitution matrix, and g=-1')
print()

sm = bs.alignment.create_substitution_matrix('ATCG', 1, -1)
_,al = bs.msa(["ATAGC", "AACC", "ATGAC"], sm, g=-1).align()

al.print()
press_enter()

print('> UPGMA')
print('     > Back to the pipeline.')
print('     > We\'ll now be running UPGMA on the resulting alignment.')

p.run_upgma()
press_enter()

print('     > Below are two different representations of the resulting phylogenetic tree.'); print()
p.text()
print()
p.phylo()
press_enter()

print('> NETWORK')
print('     > Based on the distance matrix of the alignment, we\'ll now be creating a graph/network.')
print('     > Using a cut value of 15.')
press_enter()

n = p.run_network()
print('> Some stats from the resulting network.'); print()
n.stats()
press_enter()

print('> VISUALIZATION')
print('     > Finally, we\'ll be exporting multiple files where we can visualize the previous network, with different cut values.')
print('     > Using cut values of 4, 8, 11, 21, and 40')
press_enter()

p.export_network(4)
p.export_network(8)
p.export_network(11)
p.export_network(21)
p.export_network(40)
press_enter()

print('> That\'s all, folks!')
press_enter()


