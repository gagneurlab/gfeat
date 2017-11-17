"""
Test utr3_motif_counts
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
    res = transcript.utr3_motif_counts("AC")
    assert res == 207

def test_three_letters(transcript):
    res = transcript.utr3_motif_counts("TTT")
    assert res == 151

def test_four_letters(transcript):
    res = transcript.utr3_motif_counts("ACCA")
    assert res == 14

def test_five_letters(transcript):
    res = transcript.utr3_motif_counts("GAGTT")
    assert res == 8

def test_six_letters(transcript):
    res = transcript.utr3_motif_counts("AAAAAT")
    assert res == 8

def test_seven_letters(transcript):
    res = transcript.utr3_motif_counts("TGTAAAA")
    assert res == 0
