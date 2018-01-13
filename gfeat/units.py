import pyensembl
from pysam import FastaFile
from pybedtools import BedTool
from cyvcf2 import VCF



"""Define the basic classess
"""
# basic units:
# - genes
# - transcripts
# - exons
# - introns
# - splice_junctions


# Locus:
# - contig
# - start
# - end
# - strand

class IndexUnit(object):

    @classmethod
    def from_pyensembl(cls, obj):
        raise NotImplementedError

    @classmethod
    def iter_all(cls, genome):
        raise NotImplementedError


# API independent function
#  - should be as raw as possible
def extract_sequence(interval, fasta, vcf=None):
    """function documentation - TODO
    """
    seq = fasta.fetch(str(interval.chrom), interval.start,
                       interval.stop)
    if vcf is not None:
        # TODO - modify seq according to the vcf file
        # 1. query all the positions in the vcf overlaping the interval
        # 2. For each variant:
        #    - Do a string replacement at the right places
        #      (check that the reference matches the genome sequence)
        pass
    if interval.strand == "-":
        # TODO - reverse-coplement the sequence if
        pass
    return seq

bed = BedTool("file.bed")
interval = bed[0] # get the first element in the bed file, returns an Interval object
fasta = FastaFile("file.fasta")
