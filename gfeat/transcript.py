import pyensembl
import unicodedata
# from .units import IndexUnit

# Boyer Moore string search algorithm
from pyensembl import Genome

from genome import GFGenome
from common_methods import boyer_moore_search

# class GFTranscript(IndexUnit, pyensembl.Transcript):
# Transcript class inherited from pyensembl's Transcript class
class GFTranscript(pyensembl.Transcript):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    # Counts how many codons are in the transcript's sequence
    def codon_counts(self):
        # Removing 5' UTR and 3' UTR sequences as they don't have any codons
        sequence = self.sequence.replace(self.five_prime_utr_sequence, "").replace(self.three_prime_utr_sequence, "")
        return len(sequence)/3

    # Counts how many times a given motif is presented in the 3' UTR sequence
    def utr3_motif_counts(self, pattern):
        return boyer_moore_search(self.three_prime_utr_sequence.upper(), pattern.upper())

    # Counts how many times a given motif is presented in the 5' UTR sequence
    def utr5_motif_counts(self, pattern):
        return boyer_moore_search(self.five_prime_utr_sequence.upper(), pattern.upper())

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass


    # TODO - implement
    #         - codon_counts
    #         - utr3_motif_counts()
    #            - unit-test each of them
    #               - deposit a small fasta and gtf (chr22) to the test files

# Loading testing data
# data = pyensembl.ensembl_release.EnsemblRelease(75)
# mytranscript = GFTranscript("ENST00000369985", "MYO6-001", "6", 76458926, 76629253, "+", "protein_coding", "ENSG00000196586", data)
# print(mytranscript.sequence)
# print(data.gene_ids(1, '+'))

data1 = GFGenome(reference_name='hg38_test',
                       annotation_name='hg38_chr22_test',
                       gtf_path_or_url='/Users/veronikakotova/Desktop/checking_pyensembl/gencode.v24.annotation_chr22_test.gtf',
                       transcript_fasta_paths_or_urls = '/Users/veronikakotova/Desktop/checking_pyensembl/hg38_chr22_test.fa',
                       # protein_fasta_path_or_url='/Users/veronikakotova/gfeat/tests/data/hg38_chr22.fa'
                       )
#print(dir(data1))
# data1.index()
data1.index()
# print(data1.gene_ids(22, '+'))
# print(data.transcript_ids(22, '+'))
# print(data1.gene_by_id('ENSG00000099968.17'))
# print(data1.transcript_ids_of_gene_id('ENSG00000099968.17'))
print(data1.gene_ids()[:10])
print(data1.transcript_by_id('ENST00000485631.1'))
print(data1.transcript_by_id('ENST00000485631.1').sequence)
# print(data.transcript_by_id('ENST00000415118'))
# mytranscript1 = GFTranscript('ENST00000615943.1', 'U2.14-201', '22', 10736171, 10736283, '-', 'None', 'ENSG00000277248.1', data1)
# mytranscript1 = GFTranscript("ENST00000615943.1", "U2.14-201", "ENSG00000277248.1", "U2", "None", "22:10736171-10736283")
# print(mytranscript1.sequence)

# RNA5SP46

# parse GTF and construct database of genomic features
# print(mytranscript.codon_counts())
# print(mytranscript.five_prime_utr_sequence)
# print(mytranscript.utr5_motif_counts("AC"))
