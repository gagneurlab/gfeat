"""
Test FivePrimeUTRSeq.__len__()
"""

import pytest
from tests.UTR.config_UTR import FivePrimeUTRSeq_object, FivePrimeUTRSeq_object_minus_strand


def test_len_plus(FivePrimeUTRSeq_object):
    res = len(FivePrimeUTRSeq_object)
    assert res == 8


def test_len_minus(FivePrimeUTRSeq_object_minus_strand):
    res = len(FivePrimeUTRSeq_object_minus_strand)
    assert res == 3
