"""
Test FivePrimeUTRSeq.__getitem__()
"""

import pytest
from tests.UTR.config_UTR import FivePrimeUTRSeq_object, FivePrimeUTRSeq_object_minus_strand
from pybedtools import Interval


def test_getitem_plus(FivePrimeUTRSeq_object):
    res = FivePrimeUTRSeq_object[2]
    assert res == {
        'CGCGGCCCCAAGCGGTTCCCGAGCCCAGGCCCGCGCCGAGCCCAGGTGAGCGCCCGCCCCGGGGACCCTGGTCCCCGGAACCCCGGTCCCCGCCCGCCGCCCCGCCTCTTCGCCCCGGCGCCGGGGGCCAGCGCTTGCGCTCCCAAGTCTCGGACCCCGGCCCGGCACGTTAGGGGCTGGGGGTT':
            Interval(22, 20965163, 20981358, "5' UTR", 0, "+"),
        'transcripts': ['ENST00000399163.6'],
        'exons': [
            ('CGCGGCCCCAAGCGGTTCCCGAGCCCAGGCCCGCGCCGAGCCCAG', Interval(22, 20965163, 20965207, 'ENSE00001536777.1', 0, "+")),
            (
            'GTGAGCGCCCGCCCCGGGGACCCTGGTCCCCGGAACCCCGGTCCCCGCCCGCCGCCCCGCCTCTTCGCCCCGGCGCCGGGGGCCAGCGCTTGCGCTCCCAAGTCTCGGACCCCGGCCCGGCACGTTAGGGGCTGGGGGTT',
            Interval(22, 20967805, 20967944, 'ENSE00003469659.1', 0, "+"))]
        }
