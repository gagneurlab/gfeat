"""
Test utr5_motif_counts
"""

import pytest
from tests.transcript.config import transcript


def test_two_letters(transcript):
    res = transcript.utr5_motif_counts("GG")
    assert res == 2


def test_three_letters(transcript):
    res = transcript.utr5_motif_counts("GCC")
    assert res == 1


def test_four_letters(transcript):
    res = transcript.utr5_motif_counts("GAGA")
    assert res == 0


def test_five_letters(transcript):
    res = transcript.utr5_motif_counts("ATGGC")
    assert res == 1


def test_six_letters(transcript):
    res = transcript.utr5_motif_counts("AAAAAT")
    assert res == 0
