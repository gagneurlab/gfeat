import pyensembl
import re
from itertools import product
import pandas as pd
import numpy as np


class GFTranscript(pyensembl.Transcript):
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

    def codon_counts(self):
        """
        Calculate how many codons the coding sequence has

        :return: the number of codons constituting the coding sequence
        """
        # Removing 5' UTR and 3' UTR sequences
        sequence = self.sequence.replace(self.five_prime_utr_sequence, "").replace(self.three_prime_utr_sequence, "")
        return len(sequence) / 3

    def utr3_motif_counts(self, pattern):
        """
        Calculate how many times a given motif is presented in the 3' UTR sequence

        :param pattern: string motif to be found in the 3' UTR sequence
        :type pattern: str

        :return: int, how many times a given motif is presented in the 3' UTR sequence
        """
        return len(re.findall(pattern.upper(), self.three_prime_utr_sequence.upper()))

    def utr5_motif_counts(self, pattern):
        """
        Calculate how many times a given motif is presented in the 5' UTR sequence

        :param pattern: string motif to be found in the 5' UTR sequence
        :type pattern: str

        :return: int, how many times a given motif is presented in the 5' UTR sequence
        """
        return len(re.findall(pattern.upper(), self.five_prime_utr_sequence.upper()))

    def codon_usage(self):
        """
        Calculate the frequency of all codons in the transcript

        :return: pandas.DataFrame,

                column names – 61 codons (all possible codons except for 3 stop codons),

                rows – log2 of the corresponding codon frequency
        """
        nucleobases = ['A', 'C', 'G', "T"]
        combs = [''.join(comb) for comb in product(*([nucleobases] * 3))]
        # remove stop codons
        combs.remove("TAA")
        combs.remove("TGA")
        combs.remove("TAG")

        df_codon_frequency = pd.DataFrame(np.full((1, 61), 0, dtype=int), columns=combs)

        CDS_codon_list = [self.coding_sequence.upper()[i:i + 3] for i in range(0, len(self.coding_sequence), 3)]

        for comb in combs:
            df_codon_frequency[comb][0] = CDS_codon_list.count(comb)

        # adding 1 in order to take the logarithm
        df_codon_frequency = (df_codon_frequency + 1) / len(self.coding_sequence)

        df_codon_frequency = np.log2(df_codon_frequency)

        return df_codon_frequency

    def gc_content(self, region):
        """
        Calculate the percentage of Cs and Gs in the specified region divided by 100

        :param region: 0 – coding sequence, 1 – 5'UTR sequence, 2 – 3'UTR sequence
        :type region: int

        :return: float, percentage of Cs and Gs in the specified region divided by 100
        """
        if region == 0:
            CDS = self.coding_sequence.upper()
            ratio = (CDS.count("C") + CDS.count("G")) / len(CDS)
        elif region == 1:
            CDS = self.five_prime_utr_sequence.upper()
            ratio = (CDS.count("C") + CDS.count("G")) / len(CDS)
        elif region == 2:
            CDS = self.three_prime_utr_sequence.upper()
            ratio = (CDS.count("C") + CDS.count("G")) / len(CDS)
        else:
            raise ValueError('Unknown value is provided to gc_content function.')
        return ratio

    def get_Kozak_seq(self):
        """
        Get the Kozak sequence for this transcript (6 elements upstream, start codon and 6 elements downstream)

        :return: str, Kozak sequence
        """
        utr5_seq = self.five_prime_utr_sequence.upper()
        CDS_seq = self.coding_sequence.upper()

        if len(utr5_seq) > 6:
            seq = utr5_seq[len(utr5_seq) - 6:] + CDS_seq[:9]
        else:
            seq = utr5_seq + CDS_seq[:9]
        return seq

    def get_Kozak_seq_as_df(self):
        """
        Get a line of Kozak matrix for this transcript (6 elements upstream, start codon and 6 elements downstream)

        :return: pandas.DataFrame,

                column names – first number corresponds to the base position in the Kozak sequence,
                e.g. 0A means base A 5 bases upstream from the start codon

                rows – 1 if it has the corresponding base,0 otherwise

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
        # the corresponding columns of dict_Kozak empty
        i = 15 - len(Kozak_seq)
        for c in Kozak_seq:
            dict_Kozak[str(i) + c][0] = 1
            i = i + 1

        df_Kozak = pd.DataFrame(data=dict_Kozak)

        return df_Kozak

    def get_stop_codon_context(self):
        """
        Get the stop codon context  sequence for this transcript (6 elements upstream, start codon and 6 elements downstream)

        :return: str, stop codon context  sequence
        """

        utr3_seq = self.three_prime_utr_sequence.upper()
        CDS_seq = self.coding_sequence.upper()

        seq = CDS_seq[len(CDS_seq)-9:] + utr3_seq[:6]

        return seq

    def get_stop_codon_context_as_df(self):
        """
        Get a line of stop codon context matrix for this transcript (6 elements upstream, start codon and 6 elements
        downstream)

        :return: pandas.DataFrame,

                column names – first number corresponds to the base position in the stop codon context sequence,
                e.g. 0A means base A 5 bases upstream from the start codon

                rows – 1 if it has the corresponding base, 0 otherwise
        """
        dict_stop_codon_context = {"0A": [0], "0C": [0], "0G": [0], "0T": [0],
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

        stop_codon_context = self.get_stop_codon_context()

        # apparently there are transcripts without 3'UTR regions, therefore we want to leave
        # the corresponding columns of dict_stop_codon_context empty
        i = 15 - len(stop_codon_context)
        for c in stop_codon_context:
            dict_stop_codon_context[str(i) + c][0] = 1
            i = i + 1

        df_stop_codon_context = pd.DataFrame(data=dict_stop_codon_context)

        return df_stop_codon_context

    def get_codon_pairs_frequency(self):
        """
        Calculate the frequency of 2 codons being present together in the transcript coding sequence

        :return: pandas.DataFrame,

                column names – a pair of codons, e.g. AAACAA

                rows – the frequency of the corresponding 2 codons being present together in the transcript coding sequence
        """
        nucleobases = ['A', 'C', 'G', "T"]
        combs = [''.join(comb) for comb in product(*([nucleobases] * 6))]

        iter_combs = combs.copy()

        for comb in iter_combs:
            if comb.startswith("ATG") or comb.startswith("TAA") or comb.startswith("TAG") or comb.startswith("TGA") \
                                    or comb.endswith("ATG") or comb.endswith("TAA") or comb.endswith("TAG") \
                                    or comb.endswith("TGA"):
                combs.remove(comb)

        df_codon_pairs_count = pd.DataFrame(np.full((1, 3600), 0, dtype=int), columns=combs)

        CDS_codon_list = [self.coding_sequence.upper()[i:i + 6] for i in range(0, len(self.coding_sequence), 6)]

        for comb in combs:
            df_codon_pairs_count[comb][0] = CDS_codon_list.count(comb)

        df_codon_pairs_count = (df_codon_pairs_count + 1) / len(self.coding_sequence)

        return df_codon_pairs_count
