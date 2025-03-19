from collections import Counter

def most_frequent_trinucleotide(sequence):
    trinucleotides = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
    freq = Counter(trinucleotides)
    most_common = freq.most_common(1)[0]
    return most_common

genetic_code = {
    'AUG': 'Methionine', 'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
    'UUA': 'Leucine', 'UUG': 'Leucine', 'CUU': 'Leucine', 'CUC': 'Leucine',
    'CUA': 'Leucine', 'CUG': 'Leucine', 'AUU': 'Isoleucine', 'AUC': 'Isoleucine',
    'AUA': 'Isoleucine', 'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine',
    'GUG': 'Valine', 'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine',
    'UCG': 'Serine', 'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline',
    'CCG': 'Proline', 'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine',
    'ACG': 'Threonine', 'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine',
    'GCG': 'Alanine', 'UAU': 'Tyrosine', 'UAC': 'Tyrosine', 'CAU': 'Histidine',
    'CAC': 'Histidine', 'CAA': 'Glutamine', 'CAG': 'Glutamine', 'AAU': 'Asparagine',
    'AAC': 'Asparagine', 'AAA': 'Lysine', 'AAG': 'Lysine', 'GAU': 'Aspartic Acid',
    'GAC': 'Aspartic Acid', 'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
    'UGU': 'Cysteine', 'UGC': 'Cysteine', 'UGG': 'Tryptophan', 'CGU': 'Arginine',
    'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine', 'AGU': 'Serine',
    'AGC': 'Serine', 'AGA': 'Arginine', 'AGG': 'Arginine', 'GGU': 'Glycine',
    'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine', 'UAA': 'Stop',
    'UAG': 'Stop', 'UGA': 'Stop'
}

def most_frequent_amino_acid(sequence):
    trinucleotide, _ = most_frequent_trinucleotide(sequence)
    amino_acid = genetic_code.get(trinucleotide.replace("T", "U"), "Unknown")
    return amino_acid

import matplotlib.pyplot as plt

def plot_amino_acid_frequencies(sequence):
    trinucleotides = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
    amino_acids = [genetic_code.get(tri.replace("T", "U"), "Unknown") for tri in trinucleotides]

    freq = Counter(amino_acids)
    
    plt.figure(figsize=(10, 5))
    plt.bar(freq.keys(), freq.values(), color='skyblue')
    plt.xlabel('Amino Acids')
    plt.ylabel('Frequency')
    plt.title('Amino Acid Frequency Distribution')
    plt.xticks(rotation=0)
    plt.show()

def gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

mRNA_sequence = input("输入一串RNA序列:")
print(f"The most frequent trinucleotide:{most_frequent_trinucleotide(mRNA_sequence)}")
print(f"The most frequent amino acid:{most_frequent_amino_acid(mRNA_sequence) }")
plot_amino_acid_frequencies(mRNA_sequence)
print(f"gc content: {gc_content(mRNA_sequence):.2f}%")