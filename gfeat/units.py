from cyvcf2 import VCF
from collections import OrderedDict
from itertools import combinations
from gfeat.utils import reverse_complement
import pandas as pd


class VCFMutator:
    def __init__(self, BedTool, mutations_amount, vcf=None, pos_alt_ref=None):
        """
        Constructor

        :param BedTool: True if BedTool intervals are used
        :type BedTool: bool
        :param mutations_amount: True, if the number of mutations within the interval is needed
        :type mutations_amount: bool
        :param vcf: path to the vcf.gz or file opened using cyvcf2
        :type vcf: str or an "opened" file
        :param pos_alt_ref: True if positions of mutations together with alternative and reference values
         of the nucleobase are needed
        :type pos_alt_ref: bool
        """
        self.BedTool = BedTool
        self.mutations_amount = mutations_amount
        self.vcf = vcf
        self.pos_alt_ref = pos_alt_ref

    def mutate_sequence(self, interval, fasta=None, seq_whole=None):
        """
        Function that takes an interval, fasta/seq and vcf and return all possible mutated sequences and the
        positions of mutations

        :param interval: interval object
        :type interval: pybedtools.cbedtools.Interval
        :param fasta: path to the Fasta file
        :type fasta: str
        :param seq_whole: sequence that we want to mutate
        :type seq_whole: str

        :return: depends on what was the constructor's input.

            case1: list with one tuple in it [(seq,)]

            case2: list of tuples [(DNA string, tuple with variants positions, tuple with reference values,
            tuple with alternative values), ...] and the number of all mutations within the interval
            if the number of mutations is bigger than 14, returns only the not mutated sequence since the number of
            possible sequences is too large to be processed on a laptop (as for March 2018)

            case 3: list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted

        """
        if self.vcf is None:
            # case 1: only extraction of a sequence from a fasta file is required, no vcf file is supplied
            if not fasta is None:
                if interval.strand == "-":
                    return [(reverse_complement(fasta.fetch(str(interval.chrom), interval.start,
                                                            interval.stop)),)]
                else:
                    return [(fasta.fetch(str(interval.chrom), interval.start,
                                         interval.stop),)]
            else:
                return [(seq_whole,)]
        else:
            # case 2: mutated sequences together with positions, reference and alternative values of mutations
            # are required
            if self.pos_alt_ref:
                # case 2.1: fasta file is provided
                if (not fasta == None) and (seq_whole == None):
                    seq_whole = fasta.fetch(str(interval.chrom), interval.start,
                                            interval.stop)
                # case 2.2: the path to the vcf file is provided
                if type(self.vcf) == str:
                    vcf_file = VCF(self.vcf)
                # case 2.3: the opened vcf file is provided
                else:
                    vcf_file = self.vcf
                # case 2.4: interval is of type pybedtools.cbedtools.Interval
                if not self.BedTool:
                    return _mutate_sequence_Interval_vcf_with_pos_alt_ref(interval, seq_whole, vcf_file)
                # case 2.5: interval is of type pybedtools.bedtool.BedTool
                else:
                    return _mutate_sequence_BedTool_vcf_with_pos_alt_ref(interval, seq_whole, vcf_file)
            # case 3: only mutated sequences are required
            else:
                if self.mutations_amount:
                    # case 3.1: fasta file is provided
                    if (not fasta == None) and (seq_whole == None):
                        seq_whole = fasta.fetch(str(interval.chrom), interval.start,
                                                interval.stop)
                    # case 3.2: the path to the vcf file is provided
                    if type(self.vcf) == str:
                        vcf_file = VCF(self.vcf)
                    # case 3.3: the opened vcf file is provided
                    else:
                        vcf_file = self.vcf
                    # case 3.4: interval is of type pybedtools.bedtool.BedTool
                    if self.BedTool:
                        return _mutate_sequence_BedTool_vcf(interval, seq_whole, vcf_file)
                    # case 3.5: interval is of type pybedtools.cbedtools.Interval
                    else:
                        return _mutate_sequence_Interval_vcf(interval, seq_whole, vcf_file)
                else:
                    # case 3.1: fasta file is provided
                    if (not fasta == None) and (seq_whole == None):
                        seq_whole = fasta.fetch(str(interval.chrom), interval.start,
                                                interval.stop)
                    # case 3.2: the path to the vcf file is provided
                    if type(self.vcf) == str:
                        vcf_file = VCF(self.vcf)
                    # case 3.3: the opened vcf file is provided
                    else:
                        vcf_file = self.vcf
                    # case 3.4: interval is of type pybedtools.bedtool.BedTool
                    if self.BedTool:
                        return _mutate_sequence_BedTool(interval, seq_whole, vcf_file)
                    # case 3.5: interval is of type pybedtools.cbedtools.Interval
                    else:
                        return _mutate_sequence_Interval(interval, seq_whole, vcf_file)

    def mutate_codon_context(self, intervals, seqs, column_names):
        """
        Get a table with positions of variants and their types (heterozygous or homozygous) in the given intervals for
        the given sequences

        :param intervals: intervals for which we want to know the type and positions of variants
        :type intervals: list of pybedtools.Interval objects
        :param seqs: DNA sequences corresponding to the intervals
        :type seqs: list of str
        :param column_names: prefixes to be added to the column names on order to distinguish to which interval variants
                belong to
        :type column_names: list of str
        :return: pandas.DataFrame,

            column names – corresponding prefix + the position of variant relative to the interval

            rows – 1 if the variant is heterozygous, 2 of homozygous
        """
        if self.BedTool:
            raise ValueError("Sorry! This functionality does not exist yet. Feel free to add it!")
        else:
            df_codon_context = pd.DataFrame({"name":[intervals[0].name]})
            # path to the vcf is provided
            if type(self.vcf) == str:
                vcf_file = VCF(self.vcf)
            # the opened vcf file is provided
            else:
                vcf_file = self.vcf
            i = 0
            for interval in intervals:
                vcf_interval = interval.chrom + ":" + str(interval.start) + "-" + str(interval.end - 1)
                if interval.strand == "+":
                    for variant in vcf_file(vcf_interval):
                        if variant.is_snp:
                            if (seqs[0])[variant.POS - interval.start: variant.POS - interval.start + 1] != variant.REF:
                                if (variant.POS - interval.start + 1) < 0:
                                    print("The VCF reference does not match the genome sequence")
                            else:
                                if not variant.num_het:
                                    df_codon_context[column_names[i] + str(variant.POS - interval.start)] = 2
                                else:
                                    df_codon_context[column_names[i] + str(variant.POS - interval.start)] = 1
                        else:
                            print("The variant is not a SNP")
                    i = i + 1
                else:
                    for variant in vcf_file(vcf_interval):
                        if variant.is_snp:
                            if seqs[i][variant.POS - interval.start: variant.POS - interval.start + 1] != variant.REF:
                                if (variant.POS - interval.start + 1) < 0:
                                    print("The VCF reference does not match the genome sequence")
                            else:
                                if not variant.num_het:
                                    df_codon_context[column_names[i] + str(len(seqs[i]) - (variant.POS - interval.start + 1))] = 2
                                else:
                                    df_codon_context[column_names[i] + str(len(seqs[i]) - (variant.POS - interval.start + 1))] = 1
                        else:
                            print("The variant is not a SNP")
                    i = i + 1
            return df_codon_context


def _mutate_sequence_BedTool(interval, seq, vcf_file=None):
    """
    :param interval: pybedtools.BedTool interval object
    :param seq: sequence that we want to mutate
    :param vcf: path to the vcf.gz
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
    """

    if vcf_file is not None:
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
            for variant in vcf_file(chr + ":" + str(interval.start + 1) + "-" + str(interval.end)):
                if not variant.num_het:
                    variants_hom_alt[variant.POS] = variant.ALT[0]
                    variants_hom_ref[variant.POS] = variant.REF
                else:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("Variant in position %d has 2 alternate "
                                             "non-reference alleles" % variant.POS)
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


def _mutate_sequence_Interval(interval, seq, vcf_file=None):
    """
    :param interval: pybedtools.Interval interval object
    :param seq: sequence that we want to mutate
    :param vcf: path to the vcf.gz
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
    """

    if vcf_file is not None:
        variants_hom_alt = OrderedDict()
        variants_hom_ref = OrderedDict()
        variants_het_alt = OrderedDict()
        variants_het_ref = OrderedDict()

        if interval.chrom.find("chr") == -1:
            chr = "chr" + interval.chrom
        else:
            chr = interval.chrom
        try:
            for variant in vcf_file(chr + ":" + str(interval.start) + "-" + str(interval.end)):
                if not variant.num_het:
                    variants_hom_alt[variant.POS] = variant.ALT[0]
                    variants_hom_ref[variant.POS] = variant.REF
                else:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("Variant in position %d has 2 alternate "
                                             "non-reference alleles" % variant.POS)
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


def _mutate_sequence_Interval_vcf(interval, seq, vcf_file=None):
    """
    The difference from mutate_sequence_Interval is that an already "opened" vcf file is passed to the function.

    :param interval: pybedtools.Interval interval object
    :param seq: sequence that we want to mutate
    :param vcf_file: "opened" vcf file
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
            if the number of mutations is bigger than 14, returns only the not mutated sequence since the number of
            possible sequences is too large to be processed on a laptop (as for March 2018)
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
                            raise ValueError("Variant in position %d has 2 alternate "
                                             "non-reference alleles" % variant.POS)
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


def _mutate_sequence_BedTool_vcf(interval, seq, vcf_file=None):
    """
    The difference from mutate_sequence_Interval is that an already "opened" vcf file is passed to the function.

    :param interval: pybedtools.Interval interval object
    :param seq: sequence that we want to mutate
    :param vcf_file: "opened" vcf file
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
            if the number of mutations is bigger than 14, returns only the not mutated sequence since the number of
            possible sequences is too large to be processed on a laptop (as for March 2018)
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
            for variant in vcf_file(chromosome + ":" + str(interval.start + 1) + "-" + str(interval.end)):
                if not variant.num_het:
                    variants_hom_alt[variant.POS] = variant.ALT[0]
                    variants_hom_ref[variant.POS] = variant.REF
                else:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("Variant in position %d has 2 alternate "
                                             "non-reference alleles" % variant.POS)
                        variants_het_alt[variant.POS] = variant.ALT[0]
                        variants_het_ref[variant.POS] = variant.REF

            mutated = len(variants_hom_alt) + len(variants_het_alt)

            if 0 < mutated < 14:

                tupl = ()
                output = []
                for key in variants_hom_alt:
                    if seq[key - interval.start - 1: key - interval.start +
                                                     len(variants_hom_ref[key]) - 1] != variants_hom_ref[key]:
                        raise ValueError("The VCF reference does not match the genome sequence")
                    seq = seq[:(key - interval.start) - 1] + variants_hom_alt[key] + \
                          seq[(key + len(variants_hom_alt[key]) - interval.start - 1):]
                    tupl = tupl + (key,)
                    # seq = initial sequence with all homozygous variants in it
                output.append((seq, tupl))

                for n in range(len(variants_het_alt)):
                    for combs in combinations(variants_het_alt, n + 1):
                        seq_temp = seq
                        tuple_temp = tupl
                        for key in combs:
                            if seq[key - interval.start - 1: key - interval.start +
                                                             len(variants_het_ref[key]) - 1] != variants_het_ref[key]:
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


def _mutate_sequence_Interval_vcf_with_pos_alt_ref(interval, seq, vcf_file=None):
    """
    The difference from mutate_sequence_Interval is that an already "opened" vcf file is passed to the function.

    :param interval: pybedtools.Interval interval object
    :param seq: sequence that we want to mutate
    :param vcf_file: "opened" vcf file
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
            if the number of mutations is bigger than 14, returns only the not mutated sequence since the number of
            possible sequences is too large to be processed on a laptop (as for March 2018)
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
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != "." and len(variant.ALT[0]) != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("Variant in position %d has 2 alternate "
                                             "non-reference alleles" % variant.POS)
                        variants_hom_alt[variant.POS] = variant.ALT[0]
                        variants_hom_ref[variant.POS] = variant.REF
                else:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != "." and len(variant.ALT[0]) != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("Variant in position %d has 2 alternate "
                                             "non-reference alleles" % variant.POS)
                        variants_het_alt[variant.POS] = variant.ALT[0]
                        variants_het_ref[variant.POS] = variant.REF

            mutated = len(variants_hom_alt) + len(variants_het_alt)

            if 0 <= mutated < 14:

                tupl = ()
                tupl_ref = ()
                tupl_alt = ()
                output = []
                for key in variants_hom_alt:
                    if seq[key - interval.start: key - interval.start +
                                                 len(variants_hom_ref[key])] != variants_hom_ref[key]:
                        raise ValueError("The VCF reference does not match the genome sequence")
                    seq = seq[:(key - interval.start)] + variants_hom_alt[key] + \
                          seq[(key + len(variants_hom_alt[key]) - interval.start):]
                    tupl = tupl + (key,)
                    tupl_ref = tupl_ref + (variants_hom_ref[key],)
                    tupl_alt = tupl_alt + (variants_hom_alt[key],)
                    # seq = initial sequence with all homozygous variants in it
                output.append((seq, tupl, tupl_ref, tupl_alt))

                for n in range(len(variants_het_alt)):
                    for combs in combinations(variants_het_alt, n + 1):
                        seq_temp = seq
                        tuple_temp = tupl
                        tupl_ref_temp = tupl_ref
                        tupl_alt_temp = tupl_alt
                        for key in combs:
                            if seq[key - interval.start: key - interval.start +
                                                         len(variants_het_ref[key])] != variants_het_ref[key]:
                                raise ValueError("The VCF reference does not match the genome sequence")
                            seq_temp = seq_temp[:(key - interval.start)] + variants_het_alt[key] + \
                                       seq_temp[(key + len(variants_het_alt[key]) - interval.start):]
                            tupl_ref_temp = tupl_ref_temp + (variants_het_ref[key],)
                            tupl_alt_temp = tupl_alt_temp + (variants_het_alt[key],)
                        tuple_temp = tupl + combs
                        output.append(
                            (seq_temp, tuple_temp, tupl_ref_temp, tupl_alt_temp))  # Note that the tuple is unsorted

                if interval.strand == "-":
                    temp_output = []
                    for i in range(len(output)):
                        temp_output.append((reverse_complement(output[i][0]), output[i][1], output[i][2], output[i][3]))
                    output = temp_output

            else:
                if interval.strand == "-":
                    output = [(reverse_complement(seq),)]
                else:
                    output = [(seq,)]

        except Exception as e:
            print(str(e))
            if interval.strand == "-":
                output = [(reverse_complement(seq), (), (), ())]
            else:
                output = [(seq, (), (), ())]
    else:
        if interval.strand == "-":
            output = [(reverse_complement(seq),)]
        else:
            output = [(seq,)]
    return output, mutated


def _mutate_sequence_BedTool_vcf_with_pos_alt_ref(interval, seq, vcf_file=None):
    """
    The difference from mutate_sequence_Interval is that an already "opened" vcf file is passed to the function.

    :param interval: pybedtools.Interval interval object
    :param seq: sequence that we want to mutate
    :param vcf_file: "opened" vcf file
    :return list of tuples [(DNA string, variant ids), ...] consisting of a DNA sequence and
            the positions of variants that were substituted
            if the number of mutations is bigger than 14, returns only the not mutated sequence since the number of
            possible sequences is too large to be processed on a laptop (as for March 2018)
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
            for variant in vcf_file(chromosome + ":" + str(interval.start + 1) + "-" + str(interval.end)):
                if not variant.num_het:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != "." and len(variant.ALT[0]) != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("Variant in position %d has 2 alternate "
                                             "non-reference alleles" % variant.POS)
                        variants_hom_alt[variant.POS] = variant.ALT[0]
                        variants_hom_ref[variant.POS] = variant.REF
                else:
                    # note at the moment we work only with alterations that do
                    # not introduce deletions or insertions
                    if len(variant.ALT[0]) == len(variant.REF) and variant.REF != "." and len(variant.ALT[0]) != ".":
                        # assumption: 1 ALT per 1 REF
                        if len(variant.ALT) != 1:
                            raise ValueError("Variant in position %d has 2 alternate "
                                             "non-reference alleles" % variant.POS)
                        variants_het_alt[variant.POS] = variant.ALT[0]
                        variants_het_ref[variant.POS] = variant.REF

            mutated = len(variants_hom_alt) + len(variants_het_alt)

            if 0 <= mutated < 14:

                tupl = ()
                tupl_ref = ()
                tupl_alt = ()
                output = []
                for key in variants_hom_alt:
                    if seq[key - interval.start - 1: key - interval.start +
                                                     len(variants_hom_ref[key]) - 1] != variants_hom_ref[key]:
                        raise ValueError("The VCF reference does not match the genome sequence")
                    seq = seq[:(key - interval.start) - 1] + variants_hom_alt[key] + \
                          seq[(key + len(variants_hom_alt[key]) - interval.start - 1):]
                    tupl = tupl + (key,)
                    tupl_ref = tupl_ref + (variants_hom_ref[key],)
                    tupl_alt = tupl_alt + (variants_hom_alt[key],)
                    # seq = initial sequence with all homozygous variants in it
                output.append((seq, tupl, tupl_ref, tupl_alt))

                for n in range(len(variants_het_alt)):
                    for combs in combinations(variants_het_alt, n + 1):
                        seq_temp = seq
                        tuple_temp = tupl
                        tupl_ref_temp = tupl_ref
                        tupl_alt_temp = tupl_alt
                        for key in combs:
                            if seq[key - interval.start - 1: key - interval.start +
                                                             len(variants_het_ref[key]) - 1] != variants_het_ref[key]:
                                raise ValueError("The VCF reference does not match the genome sequence")
                            seq_temp = seq_temp[:(key - interval.start - 1)] + variants_het_alt[key] + \
                                       seq_temp[(key + len(variants_het_alt[key]) - interval.start - 1):]
                            tupl_ref_temp = tupl_ref_temp + (variants_het_ref[key],)
                            tupl_alt_temp = tupl_alt_temp + (variants_het_alt[key],)
                        tuple_temp = tupl + combs
                        output.append(
                            (seq_temp, tuple_temp, tupl_ref_temp, tupl_alt_temp))  # Note that the tuple is unsorted

                if interval.strand == "-":
                    temp_output = []
                    for i in range(len(output)):
                        temp_output.append((reverse_complement(output[i][0]), output[i][1], output[i][2], output[i][3]))
                    output = temp_output

            else:
                if interval.strand == "-":
                    output = [(reverse_complement(seq),)]
                else:
                    output = [(seq,)]

        except Exception as e:
            print(str(e))
            if interval.strand == "-":
                output = [(reverse_complement(seq), (), (), ())]
            else:
                output = [(seq, (), (), ())]
    else:
        if interval.strand == "-":
            output = [(reverse_complement(seq),)]
        else:
            output = [(seq,)]
    return output, mutated
