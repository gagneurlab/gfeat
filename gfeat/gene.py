import pyensembl
import re


class GFGene(pyensembl.Gene):

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass

    def coding_sequence_codon_count(self, codon):
        """
        Calculate how many times the given codon is present in all transcript sequences corresponding to this gene
        :param codon: string, codon to be found in the sequence
        :return: how many times the given codon is present in the sequence
        """
        gene_coding_sequence = ''
        for transcript in self.transcripts:
            gene_coding_sequence += " " + transcript.coding_sequence
            return len(re.findall(codon, gene_coding_sequence))
