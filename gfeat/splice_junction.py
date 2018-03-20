import pyensembl
from gfeat.units import IndexUnit


class GFSpliceJunction(IndexUnit, pyensembl.Gene):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    # pybedtools.BedTool's interval substracts 1 from the start!!!

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass
