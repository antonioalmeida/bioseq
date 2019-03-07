import bioseq

# test = bioseq.DNASeq('AGCT')
# print(test)

# seq = bioseq.read_fasta_file('examples/dna.txt')

# {v.pretty_print() for k,v in seq.items()}

# seq2 = seq['MNPJ01000022.1 Enterocytozoon hepatopenaei strain TH1 scaffold_00014, whole genome shotgun sequence']
# #seq2.save_to_file('test.txt')

# seq3 = bioseq.DNASeq.load_from_file('test.txt')

# seq3.pretty_print()

def press_enter():
    print()
    input('Press ENTER to continue...')
    for i in range(0,5): print()

valid_dna = 'ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA'

print('Welcome to the bioseq package demo.')
print('We\'ll be using the sample DNA sequence - ' + valid_dna)
press_enter()

print('Creating and \'pretty printing\' sequence...')
dna = bioseq.DNASeq(valid_dna)
dna.pretty_print()
press_enter()

print('Computing symbol frequency...')
for k,v in dna.symbol_frequency().items():
    print('    ' + k + ' : ' + str(v) + ' occurences')
press_enter()

print('Computing GC-content...')
print('    GC-content corresponds to %.2f%% of the sequence' % (dna.gc_content()*100))
press_enter()

print('Computing transcription...')
print('The result is an instance of RNASeq.')
rna = dna.transcription()
rna.pretty_print()
press_enter()

print('Computing translation...')
print()
print('The result is the same if we start with the original DNA sequence or the result of the previous transcription execution')
aa1 = dna.translation()
aa2 = rna.translation()

print('Result of translating based on DNA sequence:') 
print('    ' + str(aa1))
print('Result of translating based on RNA sequence:')
print('    ' + str(aa2))
press_enter()
