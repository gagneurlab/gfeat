import pyensembl
import unicodedata
import re
from itertools import product
import pandas as pd
import numpy as np
# from .units import IndexUnit


# class GFTranscript(IndexUnit, pyensembl.Transcript):
# Transcript class inherited from pyensembl's Transcript class
class GFTranscript(pyensembl.Transcript):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    def CDS(self):
        """
        TODO
        :return:
        """
        # The CDS includes the start and stop codon
        sequence = self.sequence.replace(self.five_prime_utr_sequence, "").replace(self.three_prime_utr_sequence, "")
        return sequence

    def codon_counts(self):
        """

        :return: the number of codons found in the CDS
        """
        # Removing 5' UTR and 3' UTR sequences
        # Important: pyensembl 1.1.0 does not wprk correctly with user GTF and Fasta files
        sequence = self.sequence.replace(self.five_prime_utr_sequence, "").replace(self.three_prime_utr_sequence, "")
        return len(sequence)/3

    def utr3_motif_counts(self, pattern):
        """

        :param pattern: string motif to be found in the 3' UTR sequence
        :return: how many times a given motif is presented in the 3' UTR sequence
        """
        return len(re.findall(pattern.upper(), self.three_prime_utr_sequence.upper()))

    def utr5_motif_counts(self, pattern):
        """

        :param pattern: string motif to be found in the 5' UTR sequence
        :return: how many times a given motif is presented in the 5' UTR sequence
        """
        return len(re.findall(pattern.upper(), self.five_prime_utr_sequence.upper()))

    def codon_usage(self):
        """
        TODO
        :param pattern: string motif to be found in the 5' UTR sequence
        :return: how many times a given motif is presented in the 5' UTR sequence
        """

        nucleobases = ['A', 'C', 'G', "T"]
        combs = [''.join(comb) for comb in product(*([nucleobases] * 3))]
        # remove stop codons
        combs.remove("TAA")
        combs.remove("TGA")
        combs.remove("TAG")

        df_codon_frequency = pd.DataFrame(np.full((1, 61), 0, dtype=int), columns=combs)

        CDS_codon_list = [self.CDS()[i:i+3] for i in range(0, len(self.CDS()), 3)]

        # adding 1 in order to take the logarithm
        for comb in combs:
            df_codon_frequency[comb][0] = CDS_codon_list.count(comb) + 1

        df_codon_frequency = df_codon_frequency/len(self.CDS())

        df_codon_frequency = np.log2(df_codon_frequency)

        return df_codon_frequency

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass
