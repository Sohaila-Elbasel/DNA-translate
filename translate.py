# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:58:47 2020

@author: Sohaila-Elbasel
"""

def read_file(input_file):
    """ Reads and Returns a file with special characters removed """
    with open(input_file, "r") as f:
        seq = f.read()
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    
    return seq

#input_file = "dna.txt"
#file = open(input_file, "r")
#seq = file.read()
#seq = seq.replace("\n", "")
#seq = seq.replace("\r", "")
#print(seq[40:50])


def translate(seq):
    """Translate a string containing a nucleotide sequence into a string 
    containing the corresponding sequence of amino acids . Nucleotides are 
    translated in triplets using the table dictionary; each amino acid is 
    encoded with a string of length 1. """
    
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
    
    protein = ""
    if not(len(seq) % 3):
        for i in range(0, len(seq), 3):
            protein += table[seq[i : i+3]]
    return protein

dna_file = "dna.txt"
protein_file = "protein.txt"
seq = read_file(dna_file)
prt = read_file(protein_file)
print(prt == translate(seq[20:935]))