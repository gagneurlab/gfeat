import pyensembl
from .units import IndexUnit


class GFGene(IndexUnit, pyensembl.Gene):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    # Determining the possible k-mers
    # def find_kmers(string, k):
    #     kmers = []
    #     n = len(string)
    #
    #     for i in range(0, n - k + 1):
    #         kmers.append(string[i:i + k])
    #
    #     return kmers

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass


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
