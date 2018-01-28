import pyensembl
import unicodedata
import re
# from .units import IndexUnit

from pyensembl import Genome
from gfeat.genome import GFGenome


# class GFTranscript(IndexUnit, pyensembl.Transcript):
# Transcript class inherited from pyensembl's Transcript class
class GFTranscript(pyensembl.Transcript):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    def codon_counts(self):
        """

        :return: the number of codons found in the transcript's sequence
        """
        # Removing 5' UTR and 3' UTR sequences as they don't have any codons
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

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass
