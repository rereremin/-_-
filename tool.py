"""Модуль содержит классы биологических последовательностей."""


class BioSequence:
    """
    Родительский класс с общими методами для
    биологических последовательностей.
    """
    def __init__(self, sequence):
        """Конкструктор принимает последовательность."""
        self.sequence = sequence

    def get_length(self):
        """Метод для вычисления длины последовательности."""
        return len(self.sequence)

    def reverse_sequence(self):
        """Метод для вывода последовательности в обратном направлении"""
        return self.sequence[::-1]

    def count_unique_elements(self, element):
        """Метод для подсчета уникальных элементов в последовательности."""
        unique_element = sequence.count(element)
        return unique_element
    

class DNA(BioSequence):
    """Класс ДНК."""
    def __init__(self, sequence):
        """Конструктор принимает последовательность ДНК."""
        super().__init__(sequence)

    def molecular_weight(self):
        """Метод для вычисления молекулярной массы (в г/моль)."""
        dna_base_weights = {'A': 313.21, 'C': 289.18, 'G': 329.21, 'T': 304.2}
        total_weight = sum(dna_base_weights[base] for base in self.sequence)
        return total_weight

    def get_complementary_strand(self):
        """
        Метод для вывода комплиментарной цепи ДНК.
        Осуществляется по правилу Чаргаффа.
        """
        complementary_bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        complementary_strand = ''.join(complementary_bases[base] for base in self.sequence)
        return complementary_strand
    
    def get_reverse_complementary_strand(self):
        """Метод выводит последовательность, обратную комплиментарной."""
        reverse_complementary_strand = self.reverse_sequence()
        reverse_complementary_strand = DNA(reverse_complementary_strand).get_complementary_strand()
        return reverse_complementary_strand
    
    def get_mRNA_transcript(self):
        """Метод выводит мРНК, полученную на данной ДНК."""
        transcript = sequence.replace('T', 'U')
        return transcript

    def calculate_melting_temperature(self):
        """
        Метод подсчитывает температуру плавления ДНК.
        В данном методы есть верхнее ограничение в 50 пар нуклеотидов.
        Оно выбрано из соображений дизайна праймеров в ПЦР, где длинные
        фрагменты ДНК не встречаются.
        """
        if DNA(self.sequence).get_length() <= 50:
            gc_content = self.count_unique_elements('G') + self.count_unique_elements('C')
            melting_temperature = 4 * gc_content + 2 * (len(self.sequence) - gc_content)
            return melting_temperature
        return 'more than 100'

    def count_hydrogen_bonds(self):
        """Метод для подсчета количества водородных связей между цепями ДНК."""
        count = 0
        for i in range(len(self.sequence)):
            count += 3 if self.sequence[i] == 'G' or self.sequence[i] == 'C' else 2
        return count


class RNA(BioSequence):
    """Класс РНК."""
    def __init__(self, sequence):
        """Конструктор принимает последовательность РНК."""
        super().__init__(sequence)

    def molecular_weight(self):
        """Метод для вычисления молекулярной массы (в г/моль)."""
        rna_base_weights = {'A': 329.21, 'C': 305.18, 'G': 345.21, 'U': 302.16}
        total_weight = sum(rna_base_weights[base] for base in self.sequence)
        return total_weight

    def get_dna_sequence(self):
        """Получение кодирующей ДНК, с которой была образована текущая РНК."""
        dna_sequence = self.sequence.replace("U", "T")
        return dna_sequence

    def translate(self):
        """Метод осуществляется трансляцию белка на мРНК."""
        codon_table = {
            'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
            'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
            'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        protein_sequence = ""
        codons = (self.sequence[i:i+3] for i in range(0, len(self.sequence), 3))

        for codon in codons:
            if codon in codon_table:
                protein_sequence += codon_table[codon]

        return protein_sequence


class Aminoacids(BioSequence):
    """Класс аминокислотная последовательность."""
    def __init__(self, sequence):
        """Конструктор принимает последовательность аминокислот."""
        super().__init__(sequence)

    def find_proteins(self):
        """
        Метод вычленяет белки из последовательности аминокислот.
        Старт-кодон - M (метионин).
        Стоп-кодон - *.
        """
        proteins = []
        start_codon = 'M'
        stop_codon = '*'
        start_index = 0

        while True:
            start = self.sequence.find(start_codon, start_index)
            if start == -1:
                break

            end = self.sequence.find(stop_codon, start + 1)
            if end == -1:
                break

            protein = self.sequence[start:end]
            proteins.append(protein)
            start_index = end + 1

        return proteins


class Proteins(BioSequence):
    """Класс белки."""
    def __init__(self, sequence):
        """Конструктор принимает последовательность белка."""
        super().__init__(sequence)

    def percentage_hydrophobic(self):
        """Метод вычисляет содержание гидрофобных аминокислот (%)."""
        hydrophobic_aa = ['A', 'F', 'I', 'L', 'M', 'V', 'W']
        total_aa = len(self.sequence)
        count_hydrophobic = sum(aa in hydrophobic_aa for aa in self.sequence)
        percentage = (count_hydrophobic / total_aa) * 100
        return percentage

    def percentage_aromatic(self):
        """Метод вычисляет содержание ароматических аминокислот (%)."""
        aromatic_aa = ['F', 'W', 'Y']
        total_aa = len(self.sequence)
        count_aromatic = sum(aa in aromatic_aa for aa in self.sequence)
        percentage = (count_aromatic / total_aa) * 100
        return percentage

    def percentage_sulfur(self):
        """Метод вычисляет содержание серосодержащих аминокислот (%)."""
        sulfur_aa = ['C', 'M']
        total_aa = len(self.sequence)
        count_sulfur = sum(aa in sulfur_aa for aa in self.sequence)
        percentage = (count_sulfur / total_aa) * 100
        return percentage

    def molecular_weight(self):
        """Метод для вычисления молекулярной массы (в г/моль)."""
        amino_acid_weights = {
            'A': 71.09, 'R': 156.19, 'N': 114.11, 'D': 115.09, 'C': 103.15,
            'E': 129.12, 'Q': 128.14, 'G': 57.05, 'H': 137.14, 'I': 113.16,
            'L': 113.16, 'K': 128.17, 'M': 131.19, 'F': 147.18, 'P': 97.12,
            'S': 87.08, 'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13
        }

        total_weight = sum(amino_acid_weights[aa] for aa in self.sequence)
        return total_weight


def analyze_dna_sequence(sequence):
    """Вывод информации о ДНК в файл dna.txt"""
    dna = DNA(sequence)
    with open('dna.txt', mode='w') as file:
        file.write(f'Sequence: {sequence} \n')
        file.write(f'Length: {bio.get_length()} \n')
        file.write(f'Molecular weight: {round(dna.molecular_weight(), 4)} g/mole\n')
        file.write(f'Complementary strand: {dna.get_complementary_strand()} \n')
        file.write(f'Reverse complementary starnd: {dna.get_reverse_complementary_strand()} \n') 
        file.write(f'mRNA: {dna.get_mRNA_transcript()} \n') 
        file.write(f'Melting temperature is {dna.calculate_melting_temperature()}°C.\n')  
        file.write(f'Number of hydrogen bounds: {dna.count_hydrogen_bonds()}.\n')     


def analyze_rna_sequence(sequence):
    """Вывод информации о РНК в файл rna.txt"""
    rna = RNA(sequence)
    with open('rna.txt', mode='w') as file:
        file.write(f'Sequence: {sequence} \n')
        file.write(f'Length: {bio.get_length()} \n')
        file.write(f'Molecular weight: {rna.molecular_weight()} g/mole')
        file.write(f'DNA: {rna.get_dna_sequence()} \n')
        file.write(f'Aminoacids: {rna.translate()} \n')


def analyze_protein_sequence(sequence):
    """
    Вывод информации о белках в файл protein.txt.
    Перед этим осуществляется экстракция белков из 
    последовательности аминокислот.
    """
    aa = Aminoacids(sequence)
    proteins_list = aa.find_proteins()

    with open('proteins.txt', mode='w') as file:
        for protein_cur in proteins_list:
            file.write(f'Protein: {protein_cur}\n')

            protein = Proteins(protein_cur)
            file.write(f'Number of aminoacids: {protein.get_length()}\n')
            file.write(f'Molecular weight: {protein.molecular_weight()} g/mole\n')
            file.write(f'Percentage of gydrophobic aminoacids: {round(protein.percentage_hydrophobic(), 4)}%\n')
            file.write(f'Percentage of aromatic aminoacids: {round(protein.percentage_aromatic(), 4)}%\n')
            file.write(f'Percentage of sulfur aminoacids: {round(protein.percentage_sulfur(), 4)}%\n\n')


def get_properties(sequence):
    dna_bases = {'A', 'C', 'G', 'T'}
    rna_bases = {'A', 'C', 'G', 'U'}
    protein_bases = {'A', 'C', 'D', 'E', 'F', 'G', 'H',
                     'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                     'R', 'S', 'T', 'V', 'W', 'Y', '*'}
    
    if set(sequence.upper()).issubset(dna_bases):
        return analyze_dna_sequence(sequence)
    elif set(sequence.upper()).issubset(rna_bases):
        return analyze_rna_sequence(sequence)
    elif set(sequence.upper()).issubset(protein_bases):
        return analyze_protein_sequence(sequence)
    else:
        raise ValueError("Incorrect sequence!")
    

sequence = input('Input sequence: ')
bio = BioSequence(sequence.upper())
ans = get_properties(sequence.upper())