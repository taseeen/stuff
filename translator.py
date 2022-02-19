codon_to_amino_acid = {
    "UUU": "Phenylalanine",
    "UUC": "Phenylalanine",
    "UUA": "Leucine",
    "UUG": "Leucine",
    "CUU": "Leucine",
    "CUC": "Leucine",
    "CUA": "Leucine",
    "CUG": "Leucine",
    "AUU": "Isoleucine",
    "AUC": "Isoleucine",
    "AUA": "Isoleucine",
    "AUG": "Methionine",
    "GUU": "Valine",
    "GUC": "Valine",
    "GUA": "Valine",
    "GUG": "Valine",
    "UCU": "Serine",
    "UCC": "Serine",
    "UCA": "Serine",
    "UCG": "Serine",
    "CCU": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "ACU": "Threonine",
    "ACC": "Threonine",
    "ACA": "Threonine",
    "ACG": "Threonine",
    "GCU": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "UAU": "Tyrosine",
    "UAC": "Tyrosine",
    "UAA": "STOP",
    "UAG": "STOP",
    "CAU": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "AAU": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "GAU": "Aspartic Acid",
    "GAC": "Aspartic Acid",
    "GAA": "Glutamic Acid",
    "GAG": "Glutamic Acid",
    "UGU": "Cysteine",
    "UGC": "Cysteine",
    "UGA": "STOP",
    "UGG": "Tryptophan",
    "CGU": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "AGU": "Serine",
    "AGC": "Serine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GGU": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine",
}

# uracil ._.
DNA_to_RNA_complementary_base = {
    "A": "U",
    "T": "A",
    "C": "G",
    "G": "C",
}

RNA_to_DNA_complementary_base = {
    "A": "T",
    "U": "A",
    "C": "G",
    "G": "C",
}

DNA_complementary_base = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}

RNA_complementary_base = {
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C",
}

class DNA:
    def __init__(self, sequence: str):
        self.sequence = sequence[:len(sequence)-(len(sequence)%3)]
    
    def to_mRNA(self):
        mRNA_strand = ""
        for base in self.sequence:
            mRNA_strand += DNA_to_RNA_complementary_base.get(base)

        return mRNA(mRNA_strand)
    
    def to_complementary_DNA(self):
        complementary_strand = ""
        for base in self.sequence:
            complementary_strand += DNA_complementary_base.get(base)

        return DNA(complementary_strand)

class mRNA:
    def __init__(self, sequence: str):
        self.sequence = sequence[:len(sequence)-(len(sequence)%3)]
    
    def to_DNA(self):
        DNA_strand = ""
        for base in self.sequence:
            DNA_strand += RNA_to_DNA_complementary_base.get(base)

        return DNA(DNA_strand)
    
    def get_codons(self) -> list:
        return [self.sequence[i:i+3] for i in range(0, len(self.sequence), 3)]
    
    def get_anticodons(self) -> list:
        codons = self.get_codons()
        for i, codon in enumerate(codons):
            anticodon = ""
            for base in codon:
                anticodon += RNA_complementary_base.get(base)
            codons[i] = anticodon

        return codons
    
    def get_polypeptide(self) -> list:
        sequence = self.sequence[self.sequence.find("AUG"):]
        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        polypeptide = []

        for codon in codons:
            amino_acid = codon_to_amino_acid.get(codon)
            polypeptide.append(amino_acid)

            if amino_acid == "STOP":
                return polypeptide
        
        return polypeptide

DNA_input = input("Enter a DNA sequence: ").replace("\n", "").replace(" ", "").upper()

if len(DNA_input) % 3 != 0:
    print("Warning: Sequence length is not divisible by 3")

for char in DNA_input:
    if char not in ["A", "T", "C", "G"]:
        raise SyntaxError("Invalid base entered")

DNA_strand = DNA(DNA_input)

print(DNA_strand.to_mRNA().get_polypeptide())
