from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

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

def most_frequent_trinucleotide(sequence): #the first function: identify the most common trinucleotide
    stop_codons = ["UAG", "UGA", "UAA"] #three stop codons
    trinucleotides = []
    for i in range(0, len(sequence) - 2, 3): #check three base-pairs at one time
        codon = sequence[i:i+3] #put these three together
        if codon in stop_codons:
            break  #stop at the stop_codons   stop the closeest circulation
        trinucleotides.append(codon)
    freq = Counter(trinucleotides) #use Counter to statistic and return the name and number of each element in a dictionary style
    most_common = []
    for num in range((len(freq))): #to make sure that when there is more than one most frequent trinucleotide, they all can be printed
        if freq.most_common()[num][1] == freq.most_common(1)[0][1]:
            most_common.append(freq.most_common()[num][0])
    return ','.join(most_common) #use join() function to delete bracket

    #freq.most_common(1)[0][0] #use the most_common function to return a list of tuples and pick out the target tuple that has the highest statistical magnitude, and then search the tuple's first element (trinucleotide).
    #return most_common

def most_frequent_amino_acid(sequence): #the second finction: identify whether the first function's return has the corresponding most frequent amino acid
    trinucleotide = most_frequent_trinucleotide(sequence).split(",") #create a list to record
    amino_acid = []
    for item in trinucleotide: #make sure that each trinucleotide can be identified
        amino_acid.append(genetic_code[item]) #store the value on the genetic_code based on the key-trinucleotide
    return ','.join(amino_acid) #use join() function to delete bracket

def plot_amino_acid_frequencies(sequence): #the third function: draw the barchart 
    trinucleotides = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)] #check three base-pairs together
    amino_acids = [genetic_code[tri] for tri in trinucleotides] #use key of genetic_code to find the value and store.
    freq = Counter(amino_acids) #statistic and return the name and number of each element
    plt.figure(figsize=(10, 5)) #to follow the guidence
    plt.bar(freq.keys(), freq.values(), color='skyblue') #use freq.key() as the x, freq.values() as y, and use the skyblue color
    plt.xlabel('Amino Acids') #label x axis
    plt.ylabel('Frequency') #label y axis
    plt.title('Amino Acid Frequency Distribution') #name the whole picture
    plt.show()

def GC_content(sequence):
    GC_count = sequence.count('G') + sequence.count('C') #use count function to count the number of C and G
    return (GC_count / len(sequence)) * 100

def mutation(sequence):#the additional fuction:to simulate a genetic mutation and determine whether it affects the corresponding amino acid and interupts the formation of final protein
    bases=['A', 'U', 'C', 'G'] #create an array to store the four bases of RNA
    stop_codons = ["UAG", "UGA", "UAA"]
    trinucleotides = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
    for item in stop_codons:
        if item in trinucleotides:
            stop=sequence.find(item)+1
        else:
            stop=len(sequence)
    
    n=stop#to find the first base of the stop codon,where the mutation stops to be possible to happen

    position=np.random.randint(0, n) #randomly choose a position in the sequence
    current_base=sequence[position] #store the base in this position
    bases.remove(current_base)#eliminate the base in this position from the array
    new_base = np.random.choice(bases)#randomly choose a new base from the array
    mutated_sequence = sequence[:position] + new_base + sequence[position + 1:] #replace the base in the position with the new base and store the new sequence

    codon_start = position - (position % 3)#go back to the first base of the mutated codon
    #find the original codon and mutated codon
    original_codon = sequence[codon_start:codon_start+3]
    mutated_codon = mutated_sequence[codon_start:codon_start+3]
    #then using the genetic code to find the corresponding amino acid
    original_amino_acid = genetic_code[original_codon]
    mutated_amino_acid = genetic_code[mutated_codon]

    #compare original amino acid and find out the mutation type
    if mutated_amino_acid == "Stop":
        mutation_type = "nonsense mutation"
    elif original_amino_acid == mutated_amino_acid:
        mutation_type = "silent mutation"
    elif original_amino_acid != mutated_amino_acid:
        mutation_type = "missense mutation"
    
    return mutation_type, mutated_sequence

mRNA_sequence = input("Enter a string of RNA sequences(please start with the codon AUG):").strip() #eliminates any leading or trailing spaces to avoid the fault
print(f"The most frequent trinucleotide: {most_frequent_trinucleotide(mRNA_sequence)}")
print(f"The most frequent amino acid: {most_frequent_amino_acid(mRNA_sequence)}")
plot_amino_acid_frequencies(mRNA_sequence)
print(f"GC content: {GC_content(mRNA_sequence):.2f}%") #:.2f means keep two numbers after the point
mutation_type, mutated_sequence = mutation(mRNA_sequence)#call the mutation function to find out the mutation type
print(f"After the mutation, the sequence is: {mutated_sequence} \nand this mutation is: {mutation_type}") #print the result of mutation