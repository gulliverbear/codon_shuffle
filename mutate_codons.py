'''
10-13-2018
Goal is to read in genes that have already been run through 'codonshuffle'
Also read in the wt version of each gene
For each gene check the first n codons and make sure they have been mutated from WT
If they haven't then try a new codon that will mutate

In future would be good to have my own version of codon shuffle so I don't have 
to deal with all the issues using the Burge software
'''

import os
import collections
import random

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

def mutate_codon(codon, aa_to_codon, codon_to_aa):
    '''
    Randomly selects a codon that will have the most changes while 
    preserving the amino acid
    '''
    pass

def mutate_n_codons(wt_sequence, shuffled_sequence, aa_to_codon, codon_to_aa, n):
    '''
    Given wt and shuffled sequences
    Will mutate the first n codons (not counting the start ATG codon)
    returns new sequence
    '''
    new_sequence = wt_sequence[:3] # start with the ATG
    for pos in range(3, 3*n+3, 3):
        #print(f'checking position {pos}')
        wt_codon = wt_sequence[pos:pos+3]
        shuffled_codon = shuffled_sequence[pos:pos+3]
        if wt_codon == shuffled_codon:
            shuffled_codon = mutate_codon(wt_codon, aa_to_codon, codon_to_aa)
            #print(f'mutating codon {wt_codon} to {shuffled_codon}')
        #else:
            #print('no need to mutate, shuffled is different than wt')
        new_sequence += shuffled_codon
    new_sequence += shuffled_sequence[n*3+3:]
    return new_sequence

fasta_path = os.path.join('/Users','david','Documents','autoregulation','2018-09-19_shuffle_genes')
wt_fasta = os.path.abspath(os.path.join(fasta_path, 'ssu_rps_from_gary.fa'))
shuffled_fasta = os.path.abspath(os.path.join(fasta_path, '2018-09-20_ssu_rps_shuffled.fa'))

n = 10 # the # of codons at the start of each gene to mutate

wt_header_sequence = get_header_sequence(wt_fasta)
shuffled_header_sequence = get_header_sequence(shuffled_fasta)
wt_genes = set(wt_header_sequence.keys())
shuffled_genes = set(shuffled_header_sequence.keys())
# make sure wt and shuffled genes are equal
if wt_genes != shuffled_genes:
    raise SystemError('The wt and shuffled genes do not overlap')

# generate a dict to go from codon to amino acid, and vice versa
bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_to_aa = dict(zip(codons, amino_acids)) # dna: aa
aa_to_codon = collections.defaultdict(list)
for dna, amino_acid in codon_to_aa.items():
    aa_to_codon[amino_acid].append(dna) # aa: dna

# check for restriction sites in any of them
restriction_sites = ['CCTAGG', 'TTAATTAA'] # avrii, paci

#print('gene,n_shuffled_changes,n_mutated_changes,sequence_length,%_shuffle_changes,%_mutated_changed,original_sequence,shuffled_sequence,mutated_sequence')
with open('mutated_genes.fa', 'w') as f:
    for gene in shuffled_genes:
        wt_sequence = wt_header_sequence[gene]
        shuffled_sequence = shuffled_header_sequence[gene]

        has_sites = True
        while has_sites:
            mutated_sequence = mutate_n_codons(wt_sequence, shuffled_sequence,
                                              aa_to_codon, codon_to_aa, n)

            has_sites = False
            for restriction_site in restriction_sites:
                if restriction_site in mutated_sequence:
                    has_sites = True
                    print('Found restriction site so trying again')
                    break
                    
        # sanity check they are all the same length
        if len(mutated_sequence) != len(wt_sequence) and len(shuffled_sequence) != len(wt_sequence):
            raise SystemError('the sequence lengths are different')

        # count # of changes made
        num_shuffled_changes = sum([i!=j for i,j in zip(shuffled_sequence, wt_sequence)])
        num_mutated_changes = sum([i!=j for i,j in zip(mutated_sequence, wt_sequence)])

        # write out just the mutated sequences as fasta file
        f.write(f'>{gene}\n{mutated_sequence}\n')

        # make sure protein sequence wasn't changed
        wt_protein = translate_dna(wt_sequence)
        mutated_protein = translate_dna(mutated_sequence)
        if wt_protein != mutated_protein:
            raise SystemError(f'protein sequences do not match for {gene}')
            
        # for excel output
        #print(f'{gene},{num_shuffled_changes},{num_mutated_changes},{len(wt_sequence)},{num_shuffled_changes/len(wt_sequence)*100:.2f}%,{num_mutated_changes/len(wt_sequence)*100:.2f}%,{wt_sequence},{shuffled_sequence},{mutated_sequence}')


    
