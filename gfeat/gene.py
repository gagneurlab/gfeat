import pyensembl
from gfeat.common_methods import boyer_moore_search
from .genome import GFGenome


# class GFGene(IndexUnit, pyensembl.Gene):
class GFGene(pyensembl.Gene):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass

    # Counts how many times the given codob is present in the sequence
    def coding_sequence_codon_count(self, codon):
        gene_coding_sequence = ''
        for transcript in self.transcripts:
            gene_coding_sequence += transcript.coding_sequence
            return boyer_moore_search(gene_coding_sequence, codon)
#
# data1 = GFGenome(reference_name='hg38_test',
#                        annotation_name='hg38_chr22_test',
#                        gtf_path_or_url='/Users/veronikakotova/Desktop/checking_pyensembl/gencode.v24.annotation_chr22_test.gtf',
#                        transcript_fasta_paths_or_urls = '/Users/veronikakotova/Desktop/checking_pyensembl/hg38_chr22_test.fa',
#                        # protein_fasta_path_or_url='/Users/veronikakotova/gfeat/tests/data/hg38_chr22.fa'
#                        )

# print(data1.gene_by_id('ENSG00000008735.13'))
# gene = GFGene('ENSG00000008735.13', 'MAPK8IP2', '22', 50600685, 50613981, '+', None, data1)
# print(gene.coding_sequence_codon_count('AAC'))

# data = some data #for instance, data = EnsemblRelease(77)

# transcripts = transcripts(contig, strand)
# def feature_extractor(genes, feature, kmers)
    # kmers = kmers #array
    # Linear regression
    # return array of features for each gene
# def feature_extractor(transcripts, feature, kmers)
    # kmers = kmers #array
    # Linear regression
    # return array of features for each transcript
