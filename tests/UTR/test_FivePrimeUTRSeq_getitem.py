"""
Test FivePrimeUTRSeq.__getitem__()
"""

import pytest
from tests.UTR.config_UTR import FivePrimeUTRSeq_object
from pybedtools import Interval


def test_constructor(FivePrimeUTRSeq_object):
    res = FivePrimeUTRSeq_object[0]
    assert res == {'CGCGGCCCCAAGCGGTTCCCGAGCCCAGGCCCGCGCCGAGCCCAGGTGAGCGCCCGCCCCGGGGACCCTGGTCCCCGGAACCCCGGTCCCCGCCCGCCGCCCCGCCTCTTCGCCCCGGCGCCGGGGGCCAGCGCTTGCGCTCCCAAGTCTCGGACCCCGGCCCGGCACGTTAGGGGCTGGGGGTTGGCAAGCGGGGCCAGAGGCACTGGCCGGGGTCCGAAGCTCACCTCGTCCTTCTCTCGTCCCACAACGGCCCCTTCCCGGCTCGCCGCGGGGCCCCACTGTGTGCCAGTCCCCCTTCTCGGGCTGCAAGGTCTGGCGGGGGCAACGGGACGGAGCCGGCGAGGCTGCTGGCATCGGCTTCCCAGAGAGAGGGAACTGCACGTGCAGAGGCGGGGAGAGTGGACAGAGGCTTGCCCTGAGCTCCGTGGGGGAGGGGAGAGCGCGGAGCCCCCACACTCGAAGGCAGGACTGGGATCAGATGCCACCCCTAGGACGCCTGGAAGCCACTGGGAGTCCGATCTGGGCAGGAGCCCTTGGGATGAGGCAGTGGGGGCAGGTGGCGAGAGGTGGAGGAGTGGGGAGAGGACCCTGCCCAGCCTGCACCCTGGGGGCGCCAGGAGCCCTGGGAGGGTCAGGTGGAGCAGGGGCCGCTCTGTGAGATGGGAGTAATCCAGTCCCGCCTCACAGGTCCTGGGGTCCCTGGGACAGACACAGACTCCTGCTGGGTGGTTAGCACTGGGCTCTGCTAGAACTTTTCTGCAGTGGGCAGCTCATGGGCAAGGAATAGCACTCGGGCGCTCTGCACTCTGGGGGTCCCTCCTGCCCTGCCCCTGGCAAGGACCTTTGCCTTTGAGTGGGGTGACTGGGCTCCAAGCCGATCTAGCACACAGTAGGCCGAAGGCACAGGGGCCCCACCCAGTCTCCTTCCTGGCCCAGGCCCTGGGGGATTCAGGTCCCTGGGTACTCGGGAAGGGACAGCACAGCCTGGGCAGGTCCCCAGTGGCAAGACCAGGTCTGACCTTTTCTGGGATCTGTGGCCCATCAGGTGAGGCCTGAGGTCAGCCAGGAAGGAAGCGGGCAGGGCTGAACCAGGGCTGAAACCAGGACCCAGGAAGGGACAGCCTAGTGCCCAGGGAAGGGGCAGGTGGGGCAGGTGGCAGTTGAGTTTTACAGACCTTCTTTGATGACAAGTGGGGGAAAGAGCCAAGGCACGCCCCAGTGGCTGTGGGAGCCACCGCCTCCTCCCCTGTAATGGGCACCTCGAGGACCGACAGGCGCAGGGGCCAGGCACCAAGGGTCCCTCTCCAACCTGGCCTTTTGTTTCCACCCCTTGCTCTCAGGCACCCACCCACGGCCCTTCCTGGTATCTGGACCAGACCAGGTGGGTCTCAGCCCCACCTCAAGTGAGACTTCAGGAAGAACTTCCAGAGGAAACTGGAGGATGGCAGTGGGATGGACTGTGGGGAGGTCATTTGCAAACAAGGCAGCTTATGCCTGTGCCCTGGGCCTTCCAAAGTCCACCTGTAGGCTGCATCTTCCGGAAGGCCCACCCTGTAGCCTGGAAGGCACCTGGCTACAGTGATCTCTCTGCTGTGCAGGTGTCCCCTCCCCCAGGGCACCTTGCAGGGCCTTAGGGTGATATGGGGGCTTTCTGGGCCCCTGCCCACCTTAAAACTGATCAAGGTGGGAGGGTGGGGTGGGAGGCCACAGGGATGCCCCTGCTGCTTGACCTTGATGGGACCTTGACCCCTTGCTTTCCCCATATGCTGTTGGGCATGACTGCCCACAGGCCTGGTGGGCGTTTGTGAAGGGCTGGAGTCGGGGGTAGGGGGAGCTTATGGCTGCTGGGCACAGGCGTGGCAGGGCTGTGTGCCTGGCTTCCCTTTCCCCTAACCCCATCAGGCCCGGGCTGGCCGCCCACTTCTGTTCACCACAGAGCAGGGTTTCCCTGGGGACCAGTGGACTGGCCCAGCCCCACACCGCCCCCCAGCCCCTCCTGGCCACAGCTTGAGCTGGGAGGAGGGTGGAGGCTGGGCCGCCTCCTGGACCCCTGCCCTCCCTGCCCGGCCTTCCAAGGCCCGGCAGCCTCAGTCCACTGCTGGGCCTGGAACACGGAGCAGTGGCTGCCCTGCGAGGAGGTGGGTGTCCGGAGTGCTAGCCAGGGTGGGGCTGGGCTAGGGGGAGCCTGATGTCACCCCCATTGTAGGAGGGAGAAGACCAAGGATGAGAAGGGCTTTGGCTTCCCCAAGGTGACCCCCTCTTTGAGCCTAGGCATAGTGGTCTCCCGGGAGAATGGAGGGTCCCAGAGCTCTCAGGCTGTGGGATTTGGAGAGGGGAAGGGGGGCACGAGGAGGCAGTAGAAATCAGGGAGGATTTCCTGGGCAAGGCAGCCGCGGGGCTGGACGATGGGCGGCCGAGGTTCGGGCAGAGGAACTCCACTGGAGGAGTGAGGGAGCTGAGCTTGGAGGGCTTCAAAGGGGAGGTGACCTCCAGGGCCCCACGCTCTGGGATCAGGCTCACGCCATGGTTCCCAGCTCAAATGTAAGACTCTTTCCGAAGCTCCCACTGTCTCTCTCCGTCCCTCTGTGTCTCTACCTGCCCACACGCCCCGGTCCTGATGGCTCCAGTCTCCCCTGCAGGTCCTAGAGCAGCTCCAGCAGGATGGCGGCTCCAGCGTCTCTAAGGCCTGCAGGGGGTCCAGCCCCATGGGGGGCGCCCTAGGCCTCCGACAGCTCCCCATCTGTGCTCCTGCCTGCCGGCCATCCTCAGGCCACTCGCC':
                       {'transcripts': ['ENST00000399163.6'],
                        'exons': [('CGCGGCCCCAAGCGGTTCCCGAGCCCAGGCCCGCGCCGAGCCCAG', Interval(22, 20965163, 20965207, 'ENSE00001536777.1', 0)),
                        ('GTCCTAGAGCAGCTCCAGCAGGATGGCGGCTCCAGCGTCTCTAAGGCCTGCAGGGGGTCCAGCCCCATGGGGGGCGCCCTAGGCCTCCGACAGCTCCCCATCTGTGCTCCTGCCTGCCGGCCATCCTCAGGCCACTCGCC', Interval(22, 20967805, 20967944, 'ENSE00003469659.1', 0))]
                        }}
