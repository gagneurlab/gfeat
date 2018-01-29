"""
Test codon_counts
"""

from tests.transcript.config import transcript


def test_codon_counts(transcript):
    res = transcript.codon_counts()
    assert res == 2
