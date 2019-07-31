import bioseq

def press_enter():
    print()
    input('Press ENTER to continue...')
    for i in range(0,10): print()

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
print('This returns an instance of AminoacidSeq')
print()
aa1 = dna.translation()
aa2 = rna.translation()

print('Result of translating based on DNA sequence:') 
print('    ' + str(aa1))
print('Result of translating based on RNA sequence:')
print('    ' + str(aa2))
press_enter()

print('Computing all possible proteins of the resulting aminoacid sequence...')
for i in aa1.all_proteins_rf():
    print('    ' + str(i))
press_enter()

print('Computing the DNA\'s sequence codon usage of \'A\'...')
for k,v in dna.codon_usage('A').items():
    print('    ' + k + ' : ' + str(v) + ' occurences')
press_enter()

print('Computing the DNA sequence\'s reading frames...')
for i in dna.reading_frames():
    print('    ' + str(i))
press_enter()

print('Computing all putative proteins of an DNA or RNA sequence. (all open reading frames)')
for i in dna.all_open_reading_frames():
    print('    ' + str(i))
press_enter()

print('The same as above, except this time with a minimum sequence length of 6')
orf = dna.all_open_reading_frames(6)
for i in orf:
    print('    ' + str(i))
press_enter()

print('Saving the previous open reading frames (ProteinSeqs) to separates files...')
for i,o in enumerate(orf):
    o.save_to_file('protein' + str(i) + '.txt')

print('The files\' content can be later loaded using:')
print('    bioseq.DNASeq.load_from_file(filename)')
press_enter()

print('Loading a different DNA sequence, this time from a FASTA file...')
print()
print('The result is a dictionary mapping sequence description to an instance of a subclass of BioSeq.')
fasta_dic = bioseq.read_fasta_file('examples/dna.txt') 
for k,v in fasta_dic.items():
    print(k)
    v.pretty_print()
press_enter()