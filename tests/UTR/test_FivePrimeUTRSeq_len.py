"""
Test FivePrimeUTRSeq.__len__()
"""

import pytest
from tests.UTR.config_UTR import FivePrimeUTRSeq_object


def test_constructor(FivePrimeUTRSeq_object):
    res = len(FivePrimeUTRSeq_object)
    assert res == 8
