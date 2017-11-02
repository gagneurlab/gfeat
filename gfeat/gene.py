import pyensembl
from .units import IndexUnit


class GFGene(IndexUnit, pyensembl.Gene):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass
