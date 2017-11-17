"""
Test utr5_motif_counts
"""

import pytest
import pyensembl

from gfeat.transcript import GFTranscript


@pytest.fixture
def transcript():
    data = pyensembl.ensembl_release.EnsemblRelease(75)
    test_transcript = GFTranscript("ENST00000369985", "MYO6-001", "6", 76458926, 76629253, "+", "protein_coding",
                                "ENSG00000196586", data)
    return test_transcript

def test_two_letters(transcript):
    res = transcript.utr5_motif_counts("AC")
    assert res == 12

def test_three_letters(transcript):
    res = transcript.utr5_motif_counts("GGC")
    assert res == 7

def test_four_letters(transcript):
    res = transcript.utr5_motif_counts("GAGA")
    assert res == 2

def test_five_letters(transcript):
    res = transcript.utr5_motif_counts("GAGAT")
    assert res == 1

def test_six_letters(transcript):
    res = transcript.utr5_motif_counts("AAAAAT")
    assert res == 0
