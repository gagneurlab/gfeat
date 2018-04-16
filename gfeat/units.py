import pyensembl
from pysam import FastaFile
from pybedtools import BedTool
from cyvcf2 import VCF
from collections import OrderedDict
from itertools import combinations
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


def extract_sequence_BedTool(interval, fasta, vcf=None):
    """
    :param interval: pybedtools.BedTool interval object
    :param fasta: path to the Fasta file
    :param vcf: path to the vcf.gz
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
    """
    seq = fasta.fetch(str(interval.chrom), interval.start,
                       interval.stop)

    if vcf is not None:
        variants_hom_alt = OrderedDict()
        variants_hom_ref = OrderedDict()
        variants_het_alt = OrderedDict()
        variants_het_ref = OrderedDict()

        # assumption: interval.chrom looks like "chrNo"
        for variant in VCF(vcf)(interval.chrom + ":" + str(interval.start + 1) + "-" + str(interval.end)):
            if not variant.num_het:
                variants_hom_alt[variant.POS] = variant.ALT[0]
                variants_hom_ref[variant.POS] = variant.REF
            else:
                # note at the moment we work only with alterations that do
                # not introduce deletions or insertions
                if len(variant.ALT[0]) == len(variant.REF):
                    # assumption: 1 ALT per 1 REF
                    if len(variant.ALT) != 1:
                        raise ValueError("A variant in position %d has 2 alternate "
                                         "non-reference alleles" %variant.POS)
                    variants_het_alt[variant.POS] = variant.ALT[0]
                    variants_het_ref[variant.POS] = variant.REF

        tupl = ()
        output = []
        for key in variants_hom_alt:
            if seq[key - interval.start - 1: key - interval.start - 1 +
                                             len(variants_hom_ref[key])] != variants_hom_ref[key]:
                raise ValueError("The VCF reference does not match the genome sequence")
            seq = seq[:(key - interval.start - 1)] + variants_hom_alt[key] + \
                  seq[(key + len(variants_hom_alt[key]) - interval.start - 1):]
            tupl = tupl + (key,)
            # seq = initial sequence with all homozygous variants in it
        output.append((seq, tupl))

        for n in range(len(variants_het_alt)):
            for combs in combinations(variants_het_alt, n + 1):
                seq_temp = seq
                tuple_temp = tupl
                for key in combs:
                    if seq[key - interval.start - 1: key - interval.start - 1 +
                                                     len(variants_het_ref[key])] != variants_het_ref[key]:
                        raise ValueError("The VCF reference does not match the genome sequence")
                    seq_temp = seq_temp[:(key - interval.start - 1)] + variants_het_alt[key] + \
                        seq_temp[(key + len(variants_het_alt[key]) - interval.start - 1):]
                tuple_temp = tupl + combs
                output.append((seq_temp, tuple_temp))  # Note that the tuple is unsorted

        if interval.strand == "-":
            temp_output = []
            for i in range(len(output)):
                temp_output.append((reverse_complement(output[i][0]), output[i][1]))
            output = temp_output
    else:
        if interval.strand == "-":
            output = [(reverse_complement(seq),)]
        else:
            output = [(seq,)]
    return output


def mutate_sequence_BedTool(interval, seq, vcf=None):
    """
    :param interval: pybedtools.BedTool interval object
    :param seq: sequence that we want to mutate
    :param vcf: path to the vcf.gz
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
    """

    if vcf is not None:
        variants_hom_alt = OrderedDict()
        variants_hom_ref = OrderedDict()
        variants_het_alt = OrderedDict()
        variants_het_ref = OrderedDict()

        # assumption: interval.chrom looks like "chrNo"
        if (interval.chrom).find("chr") == -1:
            chr = "chr" + interval.chrom
        else:
            chr = interval.chrom
        try:
            for variant in VCF(vcf)(chr + ":" + str(interval.start + 1) + "-" + str(interval.end)):
                if not variant.num_het:
                    variants_hom_alt[variant.POS] = variant.ALT[0]
                    variants_hom_ref[variant.POS] = variant.REF
                else:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("A variant in position %d has 2 alternate "
                                             "non-reference alleles" %variant.POS)
                        variants_het_alt[variant.POS] = variant.ALT[0]
                        variants_het_ref[variant.POS] = variant.REF

            tupl = ()
            output = []
            for key in variants_hom_alt:
                if seq[key - interval.start - 1: key - interval.start - 1 +
                                                 len(variants_hom_ref[key])] != variants_hom_ref[key]:
                    raise ValueError("The VCF reference does not match the genome sequence")
                seq = seq[:(key - interval.start - 1)] + variants_hom_alt[key] + \
                      seq[(key + len(variants_hom_alt[key]) - interval.start - 1):]
                tupl = tupl + (key,)
                # seq = initial sequence with all homozygous variants in it
            output.append((seq, tupl))

            for n in range(len(variants_het_alt)):
                for combs in combinations(variants_het_alt, n + 1):
                    seq_temp = seq
                    tuple_temp = tupl
                    for key in combs:
                        if seq[key - interval.start - 1: key - interval.start - 1 +
                                                         len(variants_het_ref[key])] != variants_het_ref[key]:
                            temp = seq[key - interval.start - 1: key - interval.start - 1 +
                                                         len(variants_het_ref[key])]
                            raise ValueError("The VCF reference does not match the genome sequence")
                        seq_temp = seq_temp[:(key - interval.start - 1)] + variants_het_alt[key] + \
                            seq_temp[(key + len(variants_het_alt[key]) - interval.start - 1):]
                    tuple_temp = tupl + combs
                    output.append((seq_temp, tuple_temp))  # Note that the tuple is unsorted

            if interval.strand == "-":
                temp_output = []
                for i in range(len(output)):
                    temp_output.append((reverse_complement(output[i][0]), output[i][1]))
                output = temp_output
        except:
            if interval.strand == "-":
                output = [(reverse_complement(seq),)]
            else:
                output = [(seq,)]
    else:
        if interval.strand == "-":
            output = [(reverse_complement(seq),)]
        else:
            output = [(seq,)]
    return output


def mutate_sequence_Interval(interval, seq, vcf=None):
    """
    :param interval: pybedtools.Interval interval object
    :param seq: sequence that we want to mutate
    :param vcf: path to the vcf.gz
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
    """

    if vcf is not None:
        variants_hom_alt = OrderedDict()
        variants_hom_ref = OrderedDict()
        variants_het_alt = OrderedDict()
        variants_het_ref = OrderedDict()

        if interval.chrom.find("chr") == -1:
            chr = "chr" + interval.chrom
        else:
            chr = interval.chrom
        try:
            for variant in VCF(vcf)(chr + ":" + str(interval.start) + "-" + str(interval.end)):
                if not variant.num_het:
                    variants_hom_alt[variant.POS] = variant.ALT[0]
                    variants_hom_ref[variant.POS] = variant.REF
                else:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("A variant in position %d has 2 alternate "
                                             "non-reference alleles" %variant.POS)
                        variants_het_alt[variant.POS] = variant.ALT[0]
                        variants_het_ref[variant.POS] = variant.REF

            tupl = ()
            output = []
            for key in variants_hom_alt:
                if seq[key - interval.start: key - interval.start +
                                                 len(variants_hom_ref[key])] != variants_hom_ref[key]:
                    raise ValueError("The VCF reference does not match the genome sequence")
                seq = seq[:(key - interval.start)] + variants_hom_alt[key] + \
                      seq[(key + len(variants_hom_alt[key]) - interval.start):]
                tupl = tupl + (key,)
                # seq = initial sequence with all homozygous variants in it
            output.append((seq, tupl))

            for n in range(len(variants_het_alt)):
                for combs in combinations(variants_het_alt, n + 1):
                    seq_temp = seq
                    tuple_temp = tupl
                    for key in combs:
                        if seq[key - interval.start: key - interval.start +
                                                         len(variants_het_ref[key])] != variants_het_ref[key]:
                            raise ValueError("The VCF reference does not match the genome sequence")
                        seq_temp = seq_temp[:(key - interval.start)] + variants_het_alt[key] + \
                            seq_temp[(key + len(variants_het_alt[key]) - interval.start):]
                    tuple_temp = tupl + combs
                    output.append((seq_temp, tuple_temp))  # Note that the tuple is unsorted

            if interval.strand == "-":
                temp_output = []
                for i in range(len(output)):
                    temp_output.append((reverse_complement(output[i][0]), output[i][1]))
                output = temp_output
        except:
            if interval.strand == "-":
                output = [(reverse_complement(seq),)]
            else:
                output = [(seq,)]
    else:
        if interval.strand == "-":
            output = [(reverse_complement(seq),)]
        else:
            output = [(seq,)]
    return output


def mutate_sequence_Interval_vcf(interval, seq, vcf_file=None):
    """
    The difference from mutate_sequence_Interval is that an already "opened" vcf file is passed to the function.

    :param interval: pybedtools.Interval interval object
    :param seq: sequence that we want to mutate
    :param vcf_file: "opened" vcf file
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
            if the number of mutations is bigger than 14, returns only the not mutated sequence since the number of
            possible sequences is to large to be processed on a laptop (as for March 2018)
    :return number of mutations
    """

    if vcf_file is not None:
        variants_hom_alt = OrderedDict()
        variants_hom_ref = OrderedDict()
        variants_het_alt = OrderedDict()
        variants_het_ref = OrderedDict()

        if interval.chrom.find("chr") == -1:
            chromosome = "chr" + interval.chrom
        else:
            chromosome = interval.chrom

        mutated = 0

        try:
            for variant in vcf_file(chromosome + ":" + str(interval.start) + "-" + str(interval.end)):
                if not variant.num_het:
                    variants_hom_alt[variant.POS] = variant.ALT[0]
                    variants_hom_ref[variant.POS] = variant.REF
                else:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("A variant in position %d has 2 alternate "
                                             "non-reference alleles" %variant.POS)
                        variants_het_alt[variant.POS] = variant.ALT[0]
                        variants_het_ref[variant.POS] = variant.REF

            mutated = len(variants_hom_alt) + len(variants_het_alt)

            if 0 < mutated < 14:

                tupl = ()
                output = []
                for key in variants_hom_alt:
                    if seq[key - interval.start: key - interval.start +
                                                     len(variants_hom_ref[key])] != variants_hom_ref[key]:
                        raise ValueError("The VCF reference does not match the genome sequence")
                    seq = seq[:(key - interval.start)] + variants_hom_alt[key] + \
                          seq[(key + len(variants_hom_alt[key]) - interval.start):]
                    tupl = tupl + (key,)
                    # seq = initial sequence with all homozygous variants in it
                output.append((seq, tupl))

                for n in range(len(variants_het_alt)):
                    for combs in combinations(variants_het_alt, n + 1):
                        seq_temp = seq
                        tuple_temp = tupl
                        for key in combs:
                            if seq[key - interval.start: key - interval.start +
                                                             len(variants_het_ref[key])] != variants_het_ref[key]:
                                raise ValueError("The VCF reference does not match the genome sequence")
                            seq_temp = seq_temp[:(key - interval.start)] + variants_het_alt[key] + \
                                seq_temp[(key + len(variants_het_alt[key]) - interval.start):]
                        tuple_temp = tupl + combs
                        output.append((seq_temp, tuple_temp))  # Note that the tuple is unsorted

                if interval.strand == "-":
                    temp_output = []
                    for i in range(len(output)):
                        temp_output.append((reverse_complement(output[i][0]), output[i][1]))
                    output = temp_output

            else:
                if interval.strand == "-":
                    output = [(reverse_complement(seq),)]
                else:
                    output = [(seq,)]

        except:
            if interval.strand == "-":
                output = [(reverse_complement(seq),)]
            else:
                output = [(seq,)]
    else:
        if interval.strand == "-":
            output = [(reverse_complement(seq),)]
        else:
            output = [(seq,)]
    return output, mutated
