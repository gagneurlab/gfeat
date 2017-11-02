import pyensembl
from .units import IndexUnit


class GFTranscript(IndexUnit, pyensembl.Transcript):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

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
