import collections
import random
import os
import itertools

def translate_dna(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

def get_header_sequence(fasta_file):
    '''
    Returns dict with {header: sequence}
    '''
    header_to_seq = {}
    with open(fasta_file) as FH:
        for line in FH:
            if line.startswith('>'):
                header = line.rstrip()[1:]
                header_to_seq[header] = ''
            else:
                header_to_seq[header] += line.rstrip()
    return header_to_seq

def read_codons(sequence, codon_to_aa):
    '''
    Reads through a DNA sequence codon by codon
    Returns dic of amino acid: [list of codons]
    '''
    codons_used = collections.defaultdict(list)
    for i in range(0,len(sequence),3):
        codon = sequence[i:i+3]
        aa = codon_to_aa[codon]
        codons_used[aa].append(codon)
    return codons_used       

def shuffle_codon_lists(codons_used):
    '''
    Given a dict that has aa:[codon list]
    It will shuffle each list
    '''
    for aa in codons_used.keys():
        random.shuffle(codons_used[aa])
    return codons_used

def gene_codon_shuffle(sequence, codon_to_aa):
    '''
    to do
    '''
    pass

def shuffle_genes(header_sequence, filename, codon_to_aa):
    '''
    '''
    with open(filename, 'w') as f:
        for header, sequence in header_sequence.items():
            shuffled_sequence = gene_codon_shuffle(sequence, codon_to_aa)
            f.write(f'{header}\n{shuffled_sequence}\n')
            
            # count # of changes made
            num_shuffled_changes = sum([i!=j for i,j in zip(shuffled_sequence, sequence)])
            print(f'# of shuffled changes = {num_shuffled_changes}')

# generate a dict to go from codon to amino acid, and vice versa
bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_to_aa = dict(zip(codons, amino_acids)) # dna: aa
aa_to_codon = collections.defaultdict(list)
for dna, amino_acid in codon_to_aa.items():
    aa_to_codon[amino_acid].append(dna) # aa: dna

# read in file here as wt_fasta
wt_header_sequence = get_header_sequence(wt_fasta)
shuffle_genes(wt_header_sequence, 'shuffled_genes.fa', codon_to_aa)
