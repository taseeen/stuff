with open("translator/sequence.txt", "r") as file:
    sequence = file.read().replace("\n", "")

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
    "UAC":  "Tyrosine",
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

complementary_base = {
    "A": "U",
    "T": "A",
    "C": "G",
    "G": "C",
}

class DNA:
    def __init__(self, sequence):
        self.sequence = sequence
    
    def get_strand(self):
        return self.sequence

    def to_mRNA(self):
        mRNA = []
        for i in self.sequence:
            mRNA.append(complementary_base.get(i))
        return "".join(mRNA)

class mRNA:
    def __init__(self, sequence):
        self.sequence = sequence
    
    def get_strand(self):
        return self.sequence
    
    def get_codons(self):
        codons = []
        sequence = list(self.sequence)
        counter = 0

        while sequence:
            if counter == 2:
                counter = 0
                codons.append("".join(sequence[:3]))
                del sequence[:3]
                continue
            counter += 1

        return codons
    
    def to_polypeptide(self):
        codons = self.get_codons()

        polypeptides = []

        for codon in codons:
            amino_acid = codon_to_amino_acid.get(codon)
            polypeptides.append(amino_acid)

            if amino_acid == "STOP":
                return " - ".join(polypeptides)
        
        return " - ".join(polypeptides)

DNA_strand = DNA(sequence)
mRNA_strand = mRNA(DNA_strand.to_mRNA())
polypeptide = mRNA_strand.to_polypeptide()

print(polypeptide)