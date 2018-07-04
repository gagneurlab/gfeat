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

    def __init__(self,
                 transcript_id: object = None,
                 transcript_name: object = None,
                 contig: object = None,
                 start: object = None,
                 end: object = None,
                 strand: object = None,
                 biotype: object = None,
                 gene_id: object = None,
                 genome: object = None,
                 copy_transcript=None):
        if copy_transcript is not None:
            transcript_id = copy_transcript.transcript_id
            transcript_name = copy_transcript.transcript_name
            contig = copy_transcript.contig
            start = copy_transcript.start
            end = copy_transcript.end
            strand = copy_transcript.strand
            biotype = copy_transcript.biotype
            gene_id = copy_transcript.gene_id
            genome = copy_transcript.genome
        super(GFTranscript, self).__init__(
            transcript_id,
            transcript_name,
            contig,
            start,
            end,
            strand,
            biotype,
            gene_id,
            genome
        )

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
        return len(sequence) / 3

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

        CDS_codon_list = [self.CDS().upper()[i:i + 3] for i in range(0, len(self.CDS()), 3)]

        # adding 1 in order to take the logarithm
        for comb in combs:
            df_codon_frequency[comb][0] = CDS_codon_list.count(comb) + 1

        df_codon_frequency = df_codon_frequency / len(self.CDS())

        df_codon_frequency = np.log2(df_codon_frequency)

        return df_codon_frequency

    def GC_content(self, region):
        """
        TODO
        :param pattern: string motif to be found in the 5' UTR sequence
        :return: how many times a given motif is presented in the 5' UTR sequence
        """
        if region == 0:
            CDS = self.CDS().upper()
            ratio = (CDS.count("C") + CDS.count("G")) / len(CDS)
        elif region == 1:
            CDS = self.five_prime_utr_sequence.upper()
            ratio = (CDS.count("C") + CDS.count("G")) / len(CDS)
        elif region == 2:
            CDS = self.three_prime_utr_sequence.upper()
            ratio = (CDS.count("C") + CDS.count("G")) / len(CDS)
        else:
            raise ValueError('Unknown value is provided to GC_content function.')
        return ratio

    def get_Kozak_seq(self):
        """
        TODO
        :param pattern: string motif to be found in the 5' UTR sequence
        :return: how many times a given motif is presented in the 5' UTR sequence
        """

        utr5_seq = self.five_prime_utr_sequence.upper()
        CDS_seq = self.CDS().upper()

        seq = utr5_seq[len(utr5_seq) - 6:] + CDS_seq[:9]

        return seq

    def get_line_Kozak_matrix(self):
        """
        TODO
        :param pattern: string motif to be found in the 5' UTR sequence
        :return: how many times a given motif is presented in the 5' UTR sequence
        """
        dict_Kozak = {"0A": [0], "0C": [0], "0G": [0], "0T": [0],
                            "1A": [0], "1C": [0], "1G": [0], "1T": [0],
                            "2A": [0], "2C": [0], "2G": [0], "2T": [0],
                            "3A": [0], "3C": [0], "3G": [0], "3T": [0],
                            "4A": [0], "4C": [0], "4G": [0], "4T": [0],
                            "5A": [0], "5C": [0], "5G": [0], "5T": [0],
                            "6A": [0], "6C": [0], "6G": [0], "6T": [0],
                            "7A": [0], "7C": [0], "7G": [0], "7T": [0],
                            "8A": [0], "8C": [0], "8G": [0], "8T": [0],
                            "9A": [0], "9C": [0], "9G": [0], "9T": [0],
                            "10A": [0], "10C": [0], "10G": [0], "10T": [0],
                            "11A": [0], "11C": [0], "11G": [0], "11T": [0],
                            "12A": [0], "12C": [0], "12G": [0], "12T": [0],
                            "13A": [0], "13C": [0], "13G": [0], "13T": [0],
                            "14A": [0], "14C": [0], "14G": [0], "14T": [0]}

        Kozak_seq = self.get_Kozak_seq()

        # apparently there are transcripts without 5'UTR regions, therefore we want to leave
        # the correcponding columns of dict_Kozak empty
        i = 15 - len(Kozak_seq)
        for c in Kozak_seq:
            dict_Kozak[str(i) + c][0] = 1
            i = i + 1

        df_Kozak = pd.DataFrame(data=dict_Kozak)

        return df_Kozak

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass
