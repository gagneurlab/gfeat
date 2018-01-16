import pyensembl
from pysam import FastaFile
from pybedtools import BedTool
from cyvcf2 import VCF
from common_methods import reverse_complement



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
    """
    :param interval: pybedtools.BedTool interval object
    :param fasta: path to the Fasta file
    :param vcf: path to the vcf.gz
    :return list of tuples [(DNA sequence string, variant ids), ...]
    """
    """function documentation - TODO
    """
    seq = fasta.fetch(str(interval.chrom), interval.start,
                       interval.stop)

    print(seq)

    if vcf is not None:
        # TODO - modify seq according to the vcf file
        variants_hom = []
        variants_het = []
        # assumption: interval.chrom look like "chrNo"
        for variant in VCF(vcf).__call__(interval.chrom + ":" + str(interval.start) + "-" + str(interval.end)):
            if not variant.num_het:
                variants_hom.append((variant.ALT[0], variant.POS))
            else:
                variants_het.append((variant.ALT[0], variant.POS))
            # variant.num_het: if heterozygous equals to 1
            # variant.ALT: WHY ARE THERE 2 OF THEM??? TODO

        tuple = ()
        output = []

        for i in range(len(variants_hom)):
                seq = seq[:(variants_hom[i][1] - interval.start - 1)] + variants_hom[i][0] + \
                      seq[(variants_hom[i][1] + len(variants_hom[i][0]) - interval.start - 1):]
                tuple = tuple + (variants_hom[i][1],)
            # seq = initial sequence with all homozygous variants in it
        output.append((seq, tuple))

        for i in range(len(variants_het)):
            seq_temp = seq
            tuple_temp = tuple
            seq_temp = seq[:variants_het[i][1]] + variants_het[i][0] \
                       + seq[(variants_het[i][1] + len(variants_het[i][0])):]
            tuple_temp = tuple + (variants_het[i][1],)
            output.append((seq_temp, tuple_temp))
            for j in range(i+1,len(variants_het)):
                seq_temp = seq_temp[:variants_het[j][1]] + variants_het[j][0] \
                           + seq_temp[(variants_het[j][1]+len(variants_het[j][0])):]
                tuple_temp = tuple_temp + (variants_het[j][1],)
                output.append((seq_temp, tuple_temp))

        # 1. query all the positions in the vcf overlaping the interval
        # 2. For each variant:
        #    - Do a string replacement at the right places
        #      (check that the reference matches the genome sequence)
        pass
    if interval.strand == "-":
        # TODO - reverse-coplement the sequence if
        pass
    # return seq
    return output

gtf = BedTool("/Users/veronikakotova/gfeat/tests/data/chr22_testing.gtf")
test_interval = gtf[0]
test_fasta = FastaFile("/Users/veronikakotova/gfeat/tests/data/chr22_testing.fa")
test_vcf = "/Users/veronikakotova/gfeat/tests/data/49470G_chr22_testing.vcf.gz"

print(extract_sequence(test_interval, test_fasta, test_vcf))
# bed = BedTool("file.bed")
# interval = bed[0] # get the first element in the bed file, returns an Interval object
# fasta = FastaFile("file.fasta")
