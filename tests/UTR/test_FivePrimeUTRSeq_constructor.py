"""
Test FivePrimeUTRSeq.__init__()
"""

import pytest
from tests.UTR.config_UTR import FivePrimeUTRSeq_object, FivePrimeUTRSeq_object_minus_strand, \
    FivePrimeUTRSeq_object_minus_strand_no_reverse


def test_constructor_plus(FivePrimeUTRSeq_object):
    res = FivePrimeUTRSeq_object.seq_exons
    assert res == [
        # ' Pyensembl error',
                   'AGTGCAAACGCAGCGCCAGACACGCCCCGCCCGCGCCTCCCCGCCCCCCCGCGCCCGCGGCCCCAAGCGGTTCCCGAGCCCAGGCCCGCGCCGAGCCCAGGTGAGCGCCCGCCCCGGGGACCCTGGTCCCCGGAACCCCGGTCCCCGCCCGCCGCCCCGCCTCTTCGCCCCGGCGCCGGGGGCCAGCGCTTGCGCTCCCAAGTCTCGGACCCCGGCCCGGCACGTTAGGGGCTGGGGGTT',
                   'CGCGGCCCCAAGCGGTTCCCGAGCCCAGGCCCGCGCCGAGCCCAGGTGAGCGCCCGCCCCGGGGACCCTGGTCCCCGGAACCCCGGTCCCCGCCCGCCGCCCCGCCTCTTCGCCCCGGCGCCGGGGGCCAGCGCTTGCGCTCCCAAGTCTCGGACCCCGGCCCGGCACGTTAGGGGCTGGGGGTT',
                   'GCCCCAAGCGGTTCCCGAGCCCAGGCCCGCGCCGAGCCCAGGTGAGCGCCCGCCCCGGGGACCCTGGTCCCCGGAACCCCGGTCCCCGCCCGCCGCCCCGCCTCTTCGCCCCGGCGCCGGGGGCCAGCGCTTGCGCTCCCAAGTCTCGGACCCCGGCCCGGCACGTTAGGGGCTGGGGGTT',
                   'AAGCGGTTCCCGAGCCCAGGCCCGCGCCGAGCCCAGGTGAGCGCCCGCCCCGGGGACCCTGGTCCCCGGAACCCCGGTCCCCGCCCGCCGCCCCGCCTCTTCGCCCCGGCGCCGGGGGCCAGCGCTTGCGCTCCCAAGTCTCGGACCCCGGCCCGGCACGTTAGGGGCTGGGGGTTGGCAAGCGGGGCCAGAGGCACTGGCCGGGGTCCGAAGCTCA',
                   'CTCCTGGCCACAGCTTGAGCTGGGAGGAGGGTGGAGGCTGGGCCGCCTCCTGGACCCCTGCCCTCCCTGCCCGGCCTTCCAAGGCCCGGCAGCCTCAGTCCACTGCTGGGCCTGGAACACGGAGCAGTGGCTGCCCTGCGAGGAGGTGGGTGTCCGGAGTGCTAGCCAGGGTGGGGCTGGGCTAGGGGGAGCCTGATGTCACCCCCATTGTAGGAGGGAGAAGACCAAGGATGAGAAGGGCTTTGGCTTCCCCAAGGTGACCCCCTCTTTGAGCCTAGGCATAGT',
                   'CAGCCTCAGTCCACTGCTGGGCCTGGAACACGGAGCAGTGGCTGCCCTGCGAGGAGGTGGGTGTCCGGAGTGCTAGCCAGGGTGGGGCTGGGCTAGGGGGAGCCTGATGTCACCCCCATTGTAGGAGGGAGAAGACCAAGGATGAGAAGGGCTTTGGCTTCCCCAAGGTGACCCCCTCTTTGAGCCTAGGCATAGT',
                   'AGTCCACTGCTGGGCCTGGAACACGGAGCAGTGGCTGCCCTGCGAGGAGGTGGGTGTCCGGAGTGCTAGCCAGGGTGGGGCTGGGCTAGGGGGAGCCTGATGTCACCCCCATTGTAGGAGGGAGAAGACCAAGGATGAGAAGGGCTTTGGCTTCCCCAAGGTGACCCCCTCTTTGAGCCTAGGCATAGT',
                   'GCTCCCACTGTCTCTCTCCGTCCCTCTGTGTCTCTACCTGCCCACACGCCCCGGTCCTGATGGCTCCAGTCTCCCCTGCAGGTCCTAGAGCAGCTCCAGCAGGATGGCGGCTCCAGCGTCTCTAAGGCCTGCAGGGGGTCCAGCCCCATGGGGGGCGCCCTAGGCCTCCGACAGCTCCCCATCTGTGCTCCTGCCTGCCGGCCATCCTCAGGCCACTCGCC']


def test_constructor_minus(FivePrimeUTRSeq_object_minus_strand):
    res = FivePrimeUTRSeq_object_minus_strand.seq_exons
    assert res == [
        # ' Pyensembl error',
                   'TTTTTTTCACAAAACATTTTATTACTATAATTGATATTTTACTTATAATTTTTGGCAACATTAATAAAATAATAAATTTCACCTGAAAGAACAAGCGAGCAAAATAGACAAGGAAATTCACAAAGGGCAATAACAAAATAGTAATCTTGACATATTCAATATATTTGTAAACAAAACAAAGTGACATTGCCTTAAAAATATATAGAGGTATCAATAGAAAAACAGAGAAAGCTCCTGAAAGAACTCACTATATTAATTTTCATTTGGTTGCAGGCACAGAAATTGACTCAAGCTAGCTCAATTCTAAGACAGAAGCAAAGCCCTTGATGATGATCGCTTTCCAGCTTTTTCATGAAAACCTAGGAAATTTAAATACTTTGAAGAGGAAGAAAAGAGTGGGGAGAGAGTAAAGTGCCTTTAAGAGGAAAAGTAGAAGTTTTTTTTTTTTTTTAAAGGGAGATCTCATTGTCAGTGGCATTTTAAGAGCCTTGAAGCTCAATGAACCAAAGGAAGTGTCAAACAATTAAGTGAGGTAGCAAGCCATGCAGGCCTAGGGGAAGGGCATTCTAAGAAAGACAACAGCATGTGCAAAGTCTTTGGATTGGGAAGAATGTGATTTGCTTAAGGAATAGCAAGGCCAGTGTGTCAAAATAGTGCATCATTGGGGAAAACAGTGGAGAAGCATGACAGAGAATGAAAAAGAGAGGAAGGAGGTAGGCAGGGGCCAGATCACTGGAATCTAC',
                   'TTTTTTTTTTTTTTAAAGGGAGATCTCATTGTCAGTGGCATTTTAAGAGCCTTGAAGCTCAATGAACCAAAGGAAGTGTCAAACAATTAAGTGAGGTAGCAAGCCATGCAGGCCTAGGGGAAGGGCATTCTAAGAAAGACAACAGCATGTGCAAAGTCTTTGGATTGGGAAGAATGTGATTTGCTTAAGGAATAGCAAGGCCAGTGTGTCAAAATAGTGCATCATTGGGGAAAACAGTGGAGAAGCATGACAGAGAATGAAAAAGAGAGGAAGGAGGTAGGCAGGGGCCAGATCACTGGAATCTACAGTTGGCAATTAATGA',
                   'AAGGGAGATCTCATTGTCAGTGGCATTTTAAGAGCCTTGAAGCTCAATGAACCAAAGGAAGTGTCAAACA']


def test_constructor_minus_no_reverse(FivePrimeUTRSeq_object_minus_strand_no_reverse):
    res = FivePrimeUTRSeq_object_minus_strand_no_reverse.seq_exons
    assert res == [
        # ' Pyensembl error',
        'GTAGATTCCAGTGATCTGGCCCCTGCCTACCTCCTTCCTCTCTTTTTCATTCTCTGTCATGCTTCTCCACTGTTTTCCCCAATGATGCACTATTTTGACACACTGGCCTTGCTATTCCTTAAGCAAATCACATTCTTCCCAATCCAAAGACTTTGCACATGCTGTTGTCTTTCTTAGAATGCCCTTCCCCTAGGCCTGCATGGCTTGCTACCTCACTTAATTGTTTGACACTTCCTTTGGTTCATTGAGCTTCAAGGCTCTTAAAATGCCACTGACAATGAGATCTCCCTTTAAAAAAAAAAAAAAACTTCTACTTTTCCTCTTAAAGGCACTTTACTCTCTCCCCACTCTTTTCTTCCTCTTCAAAGTATTTAAATTTCCTAGGTTTTCATGAAAAAGCTGGAAAGCGATCATCATCAAGGGCTTTGCTTCTGTCTTAGAATTGAGCTAGCTTGAGTCAATTTCTGTGCCTGCAACCAAATGAAAATTAATATAGTGAGTTCTTTCAGGAGCTTTCTCTGTTTTTCTATTGATACCTCTATATATTTTTAAGGCAATGTCACTTTGTTTTGTTTACAAATATATTGAATATGTCAAGATTACTATTTTGTTATTGCCCTTTGTGAATTTCCTTGTCTATTTTGCTCGCTTGTTCTTTCAGGTGAAATTTATTATTTTATTAATGTTGCCAAAAATTATAAGTAAAATATCAATTATAGTAATAAAATGTTTTGTGAAAAAAA',
        'TCATTAATTGCCAACTGTAGATTCCAGTGATCTGGCCCCTGCCTACCTCCTTCCTCTCTTTTTCATTCTCTGTCATGCTTCTCCACTGTTTTCCCCAATGATGCACTATTTTGACACACTGGCCTTGCTATTCCTTAAGCAAATCACATTCTTCCCAATCCAAAGACTTTGCACATGCTGTTGTCTTTCTTAGAATGCCCTTCCCCTAGGCCTGCATGGCTTGCTACCTCACTTAATTGTTTGACACTTCCTTTGGTTCATTGAGCTTCAAGGCTCTTAAAATGCCACTGACAATGAGATCTCCCTTTAAAAAAAAAAAAAA',
        'TGTTTGACACTTCCTTTGGTTCATTGAGCTTCAAGGCTCTTAAAATGCCACTGACAATGAGATCTCCCTT']
