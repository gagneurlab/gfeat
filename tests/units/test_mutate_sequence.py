"""
Test mutate_sequence
"""

from tests.transcript.config import interval_plus, interval_minus, fasta, vcf, vcf_none, ref_check, seq
from gfeat.units import mutate_sequence


def test_extract_sequence_no_vcf(interval_plus, seq):
    res = mutate_sequence(interval_plus, seq)
    assert res == [("CCCTGAGTCATCCTTGC",)]


def test_extract_sequence_vcf_plus(interval_plus, seq, vcf):
    res = mutate_sequence(interval_plus, seq, vcf)
    assert res == [("CCCTGAGGCGTCCTTGC", (12, 14)), ("CCCCAAGGCGTCCTTGC", (12, 14, 8)),
                   ("CCCTGAAGCGTCCTTGC", (12, 14, 11)), ("CCCTGAGGCGTCTTTGC", (12, 14, 17)),
                   ("CCCCAAAGCGTCCTTGC", (12, 14, 8, 11)), ("CCCCAAGGCGTCTTTGC", (12, 14, 8, 17)),
                   ("CCCTGAAGCGTCTTTGC", (12, 14, 11, 17)), ("CCCCAAAGCGTCTTTGC", (12, 14, 8, 11, 17))]


def test_extract_sequence_vcf_minus(interval_minus, seq, vcf):
    res = mutate_sequence(interval_minus, seq, vcf)
    assert res == [("GCAAGGACGCCTCAGGG", (12, 14)), ("GCAAGGACGCCTTGGGG", (12, 14, 8)),
                   ("GCAAGGACGCTTCAGGG", (12, 14, 11)), ("GCAAAGACGCCTCAGGG", (12, 14, 17)),
                   ("GCAAGGACGCTTTGGGG", (12, 14, 8, 11)), ("GCAAAGACGCCTTGGGG", (12, 14, 8, 17)),
                   ("GCAAAGACGCTTCAGGG", (12, 14, 11, 17)), ("GCAAAGACGCTTTGGGG", (12, 14, 8, 11, 17))]
