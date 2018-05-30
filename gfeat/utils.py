def reverse_complement(dna):
    """

    :param dna: string DNA sequence
    :return: reverse-complement of a DNA sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])
