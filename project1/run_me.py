import bioseq

test = bioseq.DNASeq('AGCT')
print(test)

print(bioseq.read_fasta_file('examples/dna.txt'))